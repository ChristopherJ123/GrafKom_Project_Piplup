import { Piplup } from "./Resource/Piplup.js";
import { Prinplup } from "./Resource/Prinplup.js";
import { Empoleon } from "./Resource/Empoleon.js";
import { Environment } from "./Resource/Environment.js";
import { Geometry } from "./Resource/Geometries.js";
import { LIBS } from "./Resource/Libs.js";

class Renderer {
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        this.canvas.width = window.innerWidth;
        this.canvas.height = window.innerHeight;

        this.gl = this.canvas.getContext("webgl", { antialias: true });
        if (!this.gl) throw new Error("WebGL not supported");

        this.shader = this.createShaderProgram();

        // Create skybox shader, buffers, and texture
        this.skyboxShader = this.createSkyboxShaderProgram();
        this.skyboxBuffers = this.createSkyboxBuffers();
        this.skyboxTexture = this.loadCubemapTexture([
            'Resource/skybox/frozenrt.png', // Positive X (Right)
            'Resource/skybox/frozenlf.png', // Negative X (Left)
            'Resource/skybox/frozenup.png', // Positive Y (Top)
            'Resource/skybox/frozendn.png', // Negative Y (Bottom)
            'Resource/skybox/frozenbk.png', // Positive Z (Back)
            'Resource/skybox/frozenft.png'  // Negative Z (Front)
        ]);

        // Init objects
        this.piplup = new Piplup(this.gl, this);
        this.Prinplup = new Prinplup(this.gl, this);
        this.empoleon = new Empoleon(this.gl, this);
        this.environment = new Environment(this.gl, this);

        // Create a static matrix for Piplup's offset
        this.piplupModelMatrix = LIBS.get_I4();
        LIBS.translateX(this.piplupModelMatrix, -6.0); // Offset 6 units to the left
        LIBS.scale(this.piplupModelMatrix, 0.9); // Make it 90% of the size

        this.prinplupModelMatrix = LIBS.get_I4();
        LIBS.translateX(this.prinplupModelMatrix, -1.0);
        LIBS.translateY(this.prinplupModelMatrix, -0.18);
        LIBS.scale(this.prinplupModelMatrix, 1.2);

        // Static matrix for Empoleon
        this.empoleonModelMatrix = LIBS.get_I4();
        LIBS.translateX(this.empoleonModelMatrix, 6.0); // Offset 6 units to the right
        LIBS.translateY(this.empoleonModelMatrix, -0.4);
        LIBS.scale(this.empoleonModelMatrix, 1.9);
        
        this.viewMatrix = LIBS.get_I4();
        LIBS.translateZ(this.viewMatrix, -16);
        LIBS.translateY(this.viewMatrix, -5);
        LIBS.translateX(this.viewMatrix, 1);
        this.projMatrix = LIBS.get_projection(45, this.canvas.width / this.canvas.height, 1, 100);
        this.animationTime = 0;
        this.keysPressed = {};

        this.initInputHandlers();
        this.startRenderLoop();
    }

    createShaderProgram() {
        const gl = this.gl;
        const vsSource = `
        attribute vec3 position;
        attribute vec3 color;
        attribute vec2 texcoord;
        attribute vec3 normal;
        uniform mat4 Mmatrix, Vmatrix, Pmatrix;
        varying vec3 vColor;
        varying vec2 vTexcoord;
        varying vec3 vNormal;
        varying vec3 v_worldPosition;
        void main(void) {
            gl_Position = Pmatrix * Vmatrix * Mmatrix * vec4(position, 1.);
            vColor = color;
            vTexcoord = texcoord;
            vNormal = mat3(Mmatrix) * normal;
            v_worldPosition = (Mmatrix * vec4(position, 1.)).xyz;
        }`;
        const fsSource = `
        precision mediump float;
        varying vec3 vColor;
        varying vec2 vTexcoord;
        varying vec3 vNormal;
        varying vec3 v_worldPosition;
        uniform sampler2D sampler;
        uniform int u_useTexture;
        uniform vec3 u_lightPosition;
        uniform float u_alpha;

        void main(void) {
            vec3 normal = normalize(vNormal);
            vec3 lightDirection = normalize(u_lightPosition - v_worldPosition);
            float diff = max(dot(normal, lightDirection), 0.0);
            vec3 diffuse = diff * vec3(1.0, 1.0, 1.0);

            vec4 baseColor;
            if (u_useTexture == 1) {
                baseColor = texture2D(sampler, vTexcoord);
            } else {
                baseColor = vec4(vColor, 1.);
            }
            vec4 litColor = vec4(baseColor.rgb * (0.9 + diffuse * 0.2), baseColor.a);
            // ambient light (0.9) + diffuse light (0.2)
            gl_FragColor = vec4(baseColor.rgb * (0.9 + diffuse * 0.2), baseColor.a);
        }`;

        const vs = this.compileShader(vsSource, gl.VERTEX_SHADER);
        const fs = this.compileShader(fsSource, gl.FRAGMENT_SHADER);

        const program = gl.createProgram();
        gl.attachShader(program, vs);
        gl.attachShader(program, fs);
        gl.linkProgram(program);

        gl.useProgram(program);

        const locations = {
            position: gl.getAttribLocation(program, "position"),
            color: gl.getAttribLocation(program, "color"),
            texcoord: gl.getAttribLocation(program, "texcoord"),
            normal: gl.getAttribLocation(program, "normal"),
            Pmatrix: gl.getUniformLocation(program, "Pmatrix"),
            Vmatrix: gl.getUniformLocation(program, "Vmatrix"),
            Mmatrix: gl.getUniformLocation(program, "Mmatrix"),
            sampler: gl.getUniformLocation(program, "sampler"),
            u_useTexture: gl.getUniformLocation(program, "u_useTexture"),
            u_lightPosition: gl.getUniformLocation(program, "u_lightPosition"),
            u_alpha: gl.getUniformLocation(program, "u_alpha")
        };

        gl.enableVertexAttribArray(locations.position);
        gl.enableVertexAttribArray(locations.color);
        gl.enableVertexAttribArray(locations.texcoord);
        gl.enableVertexAttribArray(locations.normal);

        return { program, locations };
    }

    compileShader(source, type) {
        const gl = this.gl;
        const shader = gl.createShader(type);
        gl.shaderSource(shader, source);
        gl.compileShader(shader);
        if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
            throw new Error("Shader compile error: " + gl.getShaderInfoLog(shader));
        }
        return shader;
    }

    createSkyboxShaderProgram() {
        const gl = this.gl;
        const vsSource = `
        attribute vec3 position;
        varying vec3 v_texCoord;
        uniform mat4 Vmatrix;
        uniform mat4 Pmatrix;
        void main(void) {
            // Use the vertex position as the 3D texture coordinate
            v_texCoord = position;
            
            // Set the view matrix translation to 0,0,0
            // This makes the skybox "follow" the camera
            mat4 VmatrixNoTranslation = Vmatrix;
            VmatrixNoTranslation[3][0] = 0.0;
            VmatrixNoTranslation[3][1] = 0.0;
            VmatrixNoTranslation[3][2] = 0.0;
            
            gl_Position = Pmatrix * VmatrixNoTranslation * vec4(position, 1.0);
            
            // Force the skybox to be at the far clip plane (z/w = 1.0)
            // This ensures it's always drawn behind everything else.
            gl_Position = gl_Position.xyww;
        }`;
        const fsSource = `
        precision mediump float;
        varying vec3 v_texCoord;
        uniform samplerCube u_skybox;
        void main(void) {
            // Look up the color from the cubemap
            gl_FragColor = textureCube(u_skybox, normalize(v_texCoord));
        }`;

        const vs = this.compileShader(vsSource, gl.VERTEX_SHADER);
        const fs = this.compileShader(fsSource, gl.FRAGMENT_SHADER);
        const program = gl.createProgram();
        gl.attachShader(program, vs);
        gl.attachShader(program, fs);
        gl.linkProgram(program);
        gl.useProgram(program);

        const locations = {
            position: gl.getAttribLocation(program, "position"),
            Pmatrix: gl.getUniformLocation(program, "Pmatrix"),
            Vmatrix: gl.getUniformLocation(program, "Vmatrix"),
            u_skybox: gl.getUniformLocation(program, "u_skybox")
        };

        gl.enableVertexAttribArray(locations.position);
        return { program, locations };
    }

    createSkyboxBuffers() {
        const gl = this.gl;
        const skyboxGeo = Geometry.generateSkyboxCube(1.0);

        const vertexBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(skyboxGeo.vertices), gl.STATIC_DRAW);

        const facesBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, facesBuffer);
        gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(skyboxGeo.faces), gl.STATIC_DRAW);

        return {
            vertex: vertexBuffer,
            faces: facesBuffer,
            faces_length: skyboxGeo.faces.length
        };
    }

    loadCubemapTexture(urls) {
        const gl = this.gl;
        const texture = gl.createTexture();
        gl.bindTexture(gl.TEXTURE_CUBE_MAP, texture);

        const targets = [
            gl.TEXTURE_CUBE_MAP_POSITIVE_X, gl.TEXTURE_CUBE_MAP_NEGATIVE_X,
            gl.TEXTURE_CUBE_MAP_POSITIVE_Y, gl.TEXTURE_CUBE_MAP_NEGATIVE_Y,
            gl.TEXTURE_CUBE_MAP_POSITIVE_Z, gl.TEXTURE_CUBE_MAP_NEGATIVE_Z
        ];

        let imagesLoaded = 0;

        for (let i = 0; i < 6; i++) {
            const image = new Image();
            image.crossOrigin = "anonymous"; // in case you load pics from another domain
            image.onload = () => {
                gl.bindTexture(gl.TEXTURE_CUBE_MAP, texture);
                gl.texImage2D(targets[i], 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image);
                imagesLoaded++;

                if (imagesLoaded === 6) {
                    gl.generateMipmap(gl.TEXTURE_CUBE_MAP);
                    gl.texParameteri(gl.TEXTURE_CUBE_MAP, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
                }
            };
            image.src = urls[i];
        }

        gl.texParameteri(gl.TEXTURE_CUBE_MAP, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_CUBE_MAP, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

        return texture;
    }

    loadTexture(url) {
        const gl = this.gl;
        const texture = gl.createTexture();
        gl.bindTexture(gl.TEXTURE_2D, texture);

        // Placeholder pixel
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 1, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, new Uint8Array([0, 0, 255, 255]));

        const image = new Image();
        image.onload = function () {
            gl.bindTexture(gl.TEXTURE_2D, texture);
            gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image);
            gl.generateMipmap(gl.TEXTURE_2D);
        };
        image.src = url;

        return texture;
    }

    initInputHandlers() {
        let drag = false;
        let x_prev, y_prev;
        let dX_mouse = 0, dY_mouse = 0;
        let THETA = 0, PHI = 0;
        const FRICTION = 0.15;
        const KEY_ROTATION_SPEED = 0.002;

        this.canvas.onmousedown = (e) => { drag = true; x_prev = e.pageX; y_prev = e.pageY; };
        this.canvas.onmouseup = () => { drag = false; };
        this.canvas.onmouseout = () => { drag = false; };
        this.canvas.onmousemove = (e) => {
            if (!drag) return;
            dX_mouse = (e.pageX - x_prev) * 2 * Math.PI / this.canvas.width;
            dY_mouse = (e.pageY - y_prev) * 2 * Math.PI / this.canvas.height;
            x_prev = e.pageX;
            y_prev = e.pageY;
            e.preventDefault();
        };

        const keyDownHandler = (e) => {
            this.keysPressed[e.key.toLowerCase()] = true;
        };
        const keyUpHandler = (e) => {
            this.keysPressed[e.key.toLowerCase()] = false;
        };
        window.addEventListener("keydown", keyDownHandler, false);
        window.addEventListener("keyup", keyUpHandler, false);

        this.updateEnvironmentRotation = () => {
            let dX_key = 0, dY_key = 0;

            if (this.keysPressed['a']) { dX_key += KEY_ROTATION_SPEED; }
            if (this.keysPressed['d']) { dX_key -= KEY_ROTATION_SPEED; }
            if (this.keysPressed['w']) { dY_key += KEY_ROTATION_SPEED; }
            if (this.keysPressed['s']) { dY_key -= KEY_ROTATION_SPEED; }

            if (!drag) {
                dX_mouse *= (1 - FRICTION);
                dY_mouse *= (1 - FRICTION);
            }

            const totalDX = dX_mouse + dX_key;
            const totalDY = dY_mouse + dY_key;

            THETA += totalDX;
            PHI += totalDY;

            PHI = Math.max(-Math.PI / 2 + 0.1, Math.min(Math.PI / 2 - 0.1, PHI));
            
            dX_mouse = 0;
            dY_mouse = 0;

            
            const rotationMatrix = LIBS.get_I4();
            LIBS.rotateY(rotationMatrix, THETA);
            LIBS.rotateX(rotationMatrix, PHI);

            this.environment.modelMatrix = rotationMatrix;
        };
    }

    startRenderLoop() {
        const gl = this.gl;
        gl.enable(gl.DEPTH_TEST);
        gl.depthFunc(gl.LEQUAL);
        gl.clearColor(0.0, 0.0, 0.0, 0.0);
        gl.clearDepth(1.0);
        gl.enable(gl.BLEND);
        gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
        let time = 0;

        const render = (now) => {
            this.updateEnvironmentRotation();
             time += 0.02;

            const timeInSeconds = now * 0.0008;
            const bodyBreathSpeed = 1.5;
            const bodyBreathAmount = 0;
            const bodyBreathOffset  = Math.cos(timeInSeconds * bodyBreathSpeed * Math.PI) * bodyBreathAmount;
            const bodyBreathAmountScale = 0.02;
            const bodyBreathScale = 1.0 + Math.cos(timeInSeconds * bodyBreathSpeed * Math.PI) * bodyBreathAmountScale;

            const breathCycleDuration = 2 / bodyBreathSpeed;
            const timeInCycle = (timeInSeconds % breathCycleDuration);
            const breathPhase = (timeInCycle / breathCycleDuration);
            let breathAlpha = 0.0;
            let currentBreathScale = 0.1;
            const breathStartPhase = 0.0;
            const breathEndPhase = 0.5;
            const phaseDuration = breathEndPhase - breathStartPhase;

            if (breathPhase >= breathStartPhase && breathPhase <= breathEndPhase) {
                const phaseProgress = (breathPhase - breathStartPhase) / phaseDuration;
                breathAlpha = Math.sin(phaseProgress * Math.PI) * 0.4; // puncak alpha 0.4
                currentBreathScale = 0.9 + phaseProgress * 1.5; // scale up selama hembusan 0.5 to 2.0
            }

            // const eyeBreathSpeed = 1.5;
            // const eyeBreathAmount = 0.02;
            // const eyeBreathScale = 1.0 + Math.sin(timeInSeconds * eyeBreathSpeed * Math.PI) * eyeBreathAmount;

            // const diskBreathSpeed = 1.0;
            // const diskBreathAmount = 0.01;
            // const diskBreathScale = 1.0 + Math.cos(timeInSeconds * diskBreathSpeed * Math.PI) * diskBreathAmount;

            const flapSpeed = 1.0;
            const flapAmount = 0.1; // radians
            const flapAngle = Math.sin(timeInSeconds * flapSpeed * Math.PI) * flapAmount;

            const animationValues = {
                bodyTranslate: bodyBreathOffset,
                bodyScale: bodyBreathScale,
                // eye: eyeBreathScale,
                // disk: diskBreathScale,
                flapAngle: flapAngle,
                breathAlpha: breathAlpha,
                breathScale: currentBreathScale
            };

            // 1. Update all animations
            this.environment.updateAnimation();
            this.piplup.updateAnimation(timeInSeconds * 7); // ++; faster
            this.Prinplup.updateAnimation(animationValues);
            this.empoleon.updateAnimation(animationValues);


            gl.viewport(0, 0, this.canvas.width, this.canvas.height);
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

            // SKYBOX
            gl.depthMask(false);
            gl.useProgram(this.skyboxShader.program);
            // disable attributes from the main shader that are not used
            gl.disableVertexAttribArray(this.shader.locations.color);
            gl.disableVertexAttribArray(this.shader.locations.texcoord);
            gl.disableVertexAttribArray(this.shader.locations.normal);
            // position enabled (both shaders use it)
            // this.environment.modelMatrix yang holds mouse listener
            const skyboxViewMatrix = LIBS.multiply(this.viewMatrix, this.environment.modelMatrix);

            // use skyboxViewMatrix instead of this.viewMatrix
            gl.uniformMatrix4fv(this.skyboxShader.locations.Pmatrix, false, this.projMatrix);
            gl.uniformMatrix4fv(this.skyboxShader.locations.Vmatrix, false, skyboxViewMatrix);

            gl.activeTexture(gl.TEXTURE0);
            gl.bindTexture(gl.TEXTURE_CUBE_MAP, this.skyboxTexture);
            gl.uniform1i(this.skyboxShader.locations.u_skybox, 0);

            gl.bindBuffer(gl.ARRAY_BUFFER, this.skyboxBuffers.vertex);
            gl.vertexAttribPointer(this.skyboxShader.locations.position, 3, gl.FLOAT, false, 0, 0);

            gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.skyboxBuffers.faces);

            // Draw
            gl.drawElements(gl.TRIANGLES, this.skyboxBuffers.faces_length, gl.UNSIGNED_SHORT, 0);

            gl.depthMask(true); // re-enable depth writing
            gl.bindTexture(gl.TEXTURE_CUBE_MAP, null); // Unbind the cubemap from the active unit (0)

            // 2. Draw main scene
            gl.useProgram(this.shader.program); // switch back to the main shader

            // re-enable attributes for the main shader
            gl.enableVertexAttribArray(this.shader.locations.color);
            gl.enableVertexAttribArray(this.shader.locations.texcoord);
            gl.enableVertexAttribArray(this.shader.locations.normal);

            // set main scene uniforms
            gl.uniformMatrix4fv(this.shader.locations.Pmatrix, false, this.projMatrix);
            gl.uniformMatrix4fv(this.shader.locations.Vmatrix, false, this.viewMatrix);
            gl.uniform3fv(this.shader.locations.u_lightPosition, [5, 15, 10]);


            const iceParentMatrix = this.environment.getIceIslandWorldMatrix();

            const piplupParentMatrix = LIBS.multiply(this.piplupModelMatrix, iceParentMatrix);
            this.piplup.draw(this.shader, piplupParentMatrix);

            const prinplupParentMatrix = LIBS.multiply(this.prinplupModelMatrix, iceParentMatrix);
            this.Prinplup.draw(this.shader, prinplupParentMatrix, false); // false: dance. better implement with sound

            this.empoleon.updateAnimation({
                body: Math.sin(time) * 0.1,
                flapAngle: Math.sin(time * 2) * 0.2,
                tailSwing: time
            });
            const empoleonParentMatrix = LIBS.multiply(this.empoleonModelMatrix, iceParentMatrix);
            this.empoleon.draw(this.shader, empoleonParentMatrix);

            this.environment.draw(this.shader);

            requestAnimationFrame(render);
        };
        render();
    }
}

window.addEventListener('load', () => {
    new Renderer('main-canvas');
});