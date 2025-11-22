import { Piplup } from "./Resource/Piplup.js";
import { Prinplup } from "./Resource/Prinplup.js";
import { Prinplup2 } from "./Resource/Prinplup2.js";
import { Empoleon } from "./Resource/Empoleon.js";
import { Environment } from "./Resource/Environment.js";
import { Geometry } from "./Resource/Geometries.js";
import { LIBS } from "./Resource/Libs.js";
import { ModelNode } from "./Resource/ModelNode.js";

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
        this.Prinplup2 = new Prinplup2(this.gl, this);
        this.empoleon = new Empoleon(this.gl, this);
        this.environment = new Environment(this.gl, this);

        // Create a static matrix for Piplup's offset
        this.piplupModelMatrix = LIBS.get_I4();
        LIBS.translateX(this.piplupModelMatrix, -6.0); // Offset 6 units to the left
        LIBS.scale(this.piplupModelMatrix, 0.9); // Make it 90% of the size

        this.prinplupModelMatrix = LIBS.get_I4();
        LIBS.translateX(this.prinplupModelMatrix, -0.5);
        LIBS.translateY(this.prinplupModelMatrix, 6.5);
        LIBS.scale(this.prinplupModelMatrix, 0.3);

        this.prinplup2ModelMatrix = LIBS.get_I4();
        LIBS.translateX(this.prinplup2ModelMatrix, -3.0);
        LIBS.translateY(this.prinplup2ModelMatrix, -0.18);
        LIBS.scale(this.prinplup2ModelMatrix, 1.2);

        this.chatBubbleNode = this.createChatBubble(1.8, 1.5, 1.5);
        this.bubbleMatrix = LIBS.get_I4();
        LIBS.translateX(this.bubbleMatrix, 0.2); 
        LIBS.translateY(this.bubbleMatrix, 7.2);

        this.chatBubbleNode2 = this.createChatBubble(0.3, 0.3, 0.3);
        this.bubbleMatrix2 = LIBS.get_I4();
        LIBS.translateX(this.bubbleMatrix2, -1.4); 
        LIBS.translateY(this.bubbleMatrix2, 5.8);

        this.chatBubbleNode3 = this.createChatBubble(0.2, 0.2, 0.2);
        this.bubbleMatrix3 = LIBS.get_I4();
        LIBS.translateX(this.bubbleMatrix3, -1.8); 
        LIBS.translateY(this.bubbleMatrix3, 5.2);

        // Static matrix for Empoleon
        this.empoleonCenter = [6.0, -0.4, 0.0];
        this.empoleonModelMatrix = LIBS.get_I4();
        LIBS.translateX(this.empoleonModelMatrix, this.empoleonCenter[0]); // Offset 6 units to the right
        LIBS.translateY(this.empoleonModelMatrix, this.empoleonCenter[1]);
        LIBS.scale(this.empoleonModelMatrix, 1.9);
        
        this.viewMatrix = LIBS.get_I4();
        LIBS.translateZ(this.viewMatrix, -16);
        LIBS.translateY(this.viewMatrix, -5);
        LIBS.translateX(this.viewMatrix, -0.8);
        this.projMatrix = LIBS.get_projection(45, this.canvas.width / this.canvas.height, 1, 100);
        this.animationTime = 0;
        this.keysPressed = {};

        this.initInputHandlers();
        this.startRenderLoop();
    }

    createChatBubble(a,b,c) {
        const gl = this.gl;
        const whiteColor = [1.0, 1.0, 1.0];
        const bubbleGeometry = Geometry.generateSphere(a, b, c, 148, 148, whiteColor);
        const bubbleNode = new ModelNode(gl, bubbleGeometry, null);
        bubbleNode.alpha = 0.15;
        return bubbleNode;
    }
    // woi

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
            // vec4 litColor = vec4(baseColor.rgb * (0.9 + diffuse * 0.2), baseColor.a);
            // // ambient light (0.9) + diffuse light (0.2)
            // gl_FragColor = vec4(baseColor.rgb * (0.9 + diffuse * 0.2), baseColor.a);

            vec3 finalRgb = baseColor.rgb * (0.9 + diffuse * 0.2);
            gl_FragColor = vec4(finalRgb, baseColor.a * u_alpha);
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
            const bodyBreathSpeed = 0.6;
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
            this.piplup.updateAnimation(timeInSeconds * 10); // ++; faster
            this.Prinplup.updateAnimation(animationValues);
            this.Prinplup2.updateAnimation(animationValues);
            this.empoleon.updateAnimation(animationValues);

            // 1. Define circle parameters
            const circleRadius = 4.5; // How far Piplup is from Empoleon
            const circleSpeed = 1;  // How fast Piplup circles (in radians per second)
            const circleAngle = timeInSeconds * circleSpeed;

            // 2. Calculate Piplup's new X and Z position based on the circle
            const piplupX = this.empoleonCenter[0] + circleRadius * Math.cos(circleAngle);
            const piplupZ = this.empoleonCenter[2] + circleRadius * Math.sin(circleAngle);
            // Use a fixed Y level (e.g., Empoleon's Y level)
            const piplupY = this.empoleonCenter[1] + 0.4; // Adjust Y offset as needed

            // 3. Calculate rotation to make Piplup face its walking direction
            // The velocity vector is (dx/dt, dz/dt)
            // dx/dt = -circleRadius * Math.sin(circleAngle) * circleSpeed
            // dz/dt =  circleRadius * Math.cos(circleAngle) * circleSpeed
            // The angle of the tangent is atan2(dz, dx)
            const tangentAngle = Math.atan2(
                circleRadius * Math.cos(circleAngle),
                -circleRadius * Math.sin(circleAngle)
            );

            // 4. Create Piplup's new transformation matrix from scratch
            this.piplupModelMatrix = LIBS.get_I4(); // Reset to identity matrix

            // 5. Apply transformations (Order: Scale, then Rotate, then Translate)

            // a) Scale (same as you had in the constructor)
            LIBS.scale(this.piplupModelMatrix, 0.9);

            // b) Rotate Piplup to face the tangent direction
            // We subtract PI/2 because the model's "front" is +Z,
            // but the tangentAngle is 0 when moving +X.
            LIBS.rotateY(this.piplupModelMatrix, -tangentAngle + (Math.PI / 2.0));

            // c) Translate Piplup to its new position on the circle
            LIBS.translateX(this.piplupModelMatrix, piplupX);
            LIBS.translateY(this.piplupModelMatrix, piplupY);
            LIBS.translateZ(this.piplupModelMatrix, piplupZ);


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
            this.Prinplup.draw(this.shader, prinplupParentMatrix, 0); // yang kecil

            const prinplup2WorldMatrix = LIBS.multiply(this.prinplup2ModelMatrix, iceParentMatrix);
            this.Prinplup2.draw(this.shader, prinplup2WorldMatrix, 1);
            const bubbleWorldMatrix = LIBS.multiply(this.bubbleMatrix, iceParentMatrix);
            // LIBS.translateX(bubbleWorldMatrix, -3);
            // LIBS.translateY(bubbleWorldMatrix, 5);
            // LIBS.translateZ(bubbleWorldMatrix, 0);
            this.chatBubbleNode.updateWorldMatrix(bubbleWorldMatrix);

            const bubbleWorldMatrix2 = LIBS.multiply(this.bubbleMatrix2, iceParentMatrix);
            this.chatBubbleNode2.updateWorldMatrix(bubbleWorldMatrix2);

            const bubbleWorldMatrix3 = LIBS.multiply(this.bubbleMatrix3, iceParentMatrix);
            this.chatBubbleNode3.updateWorldMatrix(bubbleWorldMatrix3);

            // gl.depthMask(false);
            this.chatBubbleNode.draw(this.shader);
            this.chatBubbleNode2.draw(this.shader);
            this.chatBubbleNode3.draw(this.shader);
            gl.depthMask(true);
            // woi
            

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