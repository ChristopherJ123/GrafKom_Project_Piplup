import { Piplup } from "./Piplup.js";
import { Prinplup } from "./Prinplup.js";
import { Empoleon } from "./Empoleon.js";
import { Environment } from "./Environment.js";
import { Geometry } from "./Geometries.js";
import { LIBS } from "./Libs.js";

class Renderer {
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        this.canvas.width = window.innerWidth;
        this.canvas.height = window.innerHeight;

        this.gl = this.canvas.getContext("webgl", { antialias: true });
        if (!this.gl) throw new Error("WebGL not supported");

        this.shader = this.createShaderProgram();

        // NEW: Create skybox shader, buffers, and texture
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

        // Piplup
        this.piplup = new Piplup(this.gl, this);

        this.Prinplup = new Prinplup(this.gl, this);

        // Empoleon
        this.empoleon = new Empoleon(this.gl, this);

        this.environment = new Environment(this.gl, this);

        // NEW: Create a static matrix for Piplup's offset
        this.piplupModelMatrix = LIBS.get_I4();
        LIBS.translateX(this.piplupModelMatrix, -6.0); // Position 5 units to the left
        LIBS.scale(this.piplupModelMatrix, 0.9); // Make it 90% of the size

        this.prinplupModelMatrix = LIBS.get_I4();
        LIBS.translateX(this.prinplupModelMatrix, -1.0);
        LIBS.translateY(this.prinplupModelMatrix, -0.04);
        LIBS.scale(this.prinplupModelMatrix, 1.2);

        // Static matrix for Empoleon
        this.empoleonModelMatrix = LIBS.get_I4();
        LIBS.translateX(this.empoleonModelMatrix, 6.0); // Position 5 units to the right
        LIBS.translateY(this.empoleonModelMatrix, -0.4); // Position 5 units to the right
        LIBS.scale(this.empoleonModelMatrix, 1.9); // Keep it at 100% size
        
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
            // Apply lighting
            vec4 litColor = vec4(baseColor.rgb * (0.9 + diffuse * 0.2), baseColor.a);
            // Apply ambient light (0.9) + diffuse light (0.2)
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

    /**
     * Creates a new shader program specifically for the skybox.
     * This shader is simple: it's not affected by lighting.
     */
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

    /**
     * Creates the vertex and face buffers for the skybox.
     */
    createSkyboxBuffers() {
        const gl = this.gl;
        const skyboxGeo = Geometry.generateSkyboxCube(1.0); // Size doesn't matter, we scale it

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

    /**
     * Loads 6 images into a single CUBE_MAP texture.
     */
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
            image.crossOrigin = "anonymous"; // In case you load from another domain
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
            // Calculate mouse delta only
            dX_mouse = (e.pageX - x_prev) * 2 * Math.PI / this.canvas.width;
            dY_mouse = (e.pageY - y_prev) * 2 * Math.PI / this.canvas.height;
            // Don't directly add to THETA/PHI here, let updateEnvironmentRotation handle it
            x_prev = e.pageX;
            y_prev = e.pageY;
            e.preventDefault();
        };

        // --- Keyboard Event Listeners (Keep as before, just track keys) ---
        const keyDownHandler = (e) => {
            this.keysPressed[e.key.toLowerCase()] = true;
        };
        const keyUpHandler = (e) => {
            this.keysPressed[e.key.toLowerCase()] = false;
        };
        window.addEventListener("keydown", keyDownHandler, false);
        window.addEventListener("keyup", keyUpHandler, false);

        // --- HAPUS FUNGSI LAMA INI ---
        // this.updateMouseRotation = () => { ... };
        // this.updateCameraRotation = () => { ... };

        // --- TAMBAHKAN FUNGSI BARU INI (Gabungan) ---
        this.updateEnvironmentRotation = () => {
            let dX_key = 0, dY_key = 0; // Keyboard delta

            // Calculate keyboard delta
            if (this.keysPressed['a']) { dX_key += KEY_ROTATION_SPEED; }
            if (this.keysPressed['d']) { dX_key -= KEY_ROTATION_SPEED; }
            if (this.keysPressed['w']) { dY_key += KEY_ROTATION_SPEED; }
            if (this.keysPressed['s']) { dY_key -= KEY_ROTATION_SPEED; }

            // Apply friction if mouse is not dragging
            if (!drag) {
                dX_mouse *= (1 - FRICTION);
                dY_mouse *= (1 - FRICTION);
            }

            // Combine deltas from mouse and keyboard
            const totalDX = dX_mouse + dX_key;
            const totalDY = dY_mouse + dY_key;

            // Update angles
            THETA += totalDX;
            PHI += totalDY;

            // Limit vertical rotation (optional but recommended)
            PHI = Math.max(-Math.PI / 2 + 0.1, Math.min(Math.PI / 2 - 0.1, PHI));

            // Reset mouse delta after applying (keyboard delta is recalculated each frame)
            dX_mouse = 0;
            dY_mouse = 0;

            // Create the final rotation matrix
            const rotationMatrix = LIBS.get_I4();
            LIBS.rotateY(rotationMatrix, THETA); // Horizontal rotation
            LIBS.rotateX(rotationMatrix, PHI);   // Vertical rotation

            // Apply the combined rotation to the environment's model matrix
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

        const render = (now) => {
            this.updateEnvironmentRotation();

            // --- CALCULATE ANIMATION VALUES ---
            const timeInSeconds = now * 0.0008;
            const bodyBreathSpeed = 1.5;
            const bodyBreathAmount = 0.02;
            const bodyBreathOffset  = Math.cos(timeInSeconds * bodyBreathSpeed * Math.PI) * bodyBreathAmount;

            // Body breathe - Scale
            const bodyBreathAmountScale = 0.02;
            const bodyBreathScale = 1.0 + Math.cos(timeInSeconds * bodyBreathSpeed * Math.PI) * bodyBreathAmountScale;

            const breathCycleDuration = 2 / bodyBreathSpeed; // Time for one inhale/exhale
            const timeInCycle = (timeInSeconds % breathCycleDuration);
            const breathPhase = (timeInCycle / breathCycleDuration); // 0 to 1
            let breathAlpha = 0.0;
            let currentBreathScale = 0.1; // Start small
            const breathStartPhase = 0.0; // Mulai hembusan saat tubuh mulai mengecil
            const breathEndPhase = 0.5;   // Akhiri hembusan saat tubuh paling kecil
            const phaseDuration = breathEndPhase - breathStartPhase; // Durasi = 0.5

            if (breathPhase >= breathStartPhase && breathPhase <= breathEndPhase) {
                const phaseProgress = (breathPhase - breathStartPhase) / phaseDuration; // 0 sampai 1
                // Gunakan sin untuk fade in/out yang mulus
                breathAlpha = Math.sin(phaseProgress * Math.PI) * 0.4; // Puncak alpha 0.4
                // Skala membesar selama hembusan
                currentBreathScale = 0.9 + phaseProgress * 1.5; // Skala dari 0.5 sampai 2.0
            }

            const eyeBreathSpeed = 1.5;
            const eyeBreathAmount = 0.02;
            const eyeBreathScale = 1.0 + Math.sin(timeInSeconds * eyeBreathSpeed * Math.PI) * eyeBreathAmount;

            const diskBreathSpeed = 1.0;
            const diskBreathAmount = 0.01;
            const diskBreathScale = 1.0 + Math.cos(timeInSeconds * diskBreathSpeed * Math.PI) * diskBreathAmount;

            const flapSpeed = 1.0;
            const flapAmount = 0.1; // Radians
            const flapAngle = Math.sin(timeInSeconds * flapSpeed * Math.PI) * flapAmount;

            const animationValues = {
                bodyTranslate: bodyBreathOffset,
                bodyScale: bodyBreathScale,
                eye: eyeBreathScale,
                disk: diskBreathScale,
                flapAngle: flapAngle,
                breathAlpha: breathAlpha,
                breathScale: currentBreathScale
            };
            // --- END ANIMATION VALUES ---

            // --- UPDATE MODELS ---
            this.environment.updateAnimation(); // Updates ice island float

            // NEW: Update Piplup's animation
            // Note: Piplup's update function just takes 'time', not 'animationValues'
            this.piplup.updateAnimation(timeInSeconds * 7); // *1.5 to make it run a bit faster

            this.Prinplup.updateAnimation(animationValues); // Updates body/hand anims

            // Empoleon update animation
            this.empoleon.updateAnimation(animationValues);


            // --- RENDER ---
            gl.viewport(0, 0, this.canvas.width, this.canvas.height);
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

            // === 1. DRAW SKYBOX ===
            gl.depthMask(false); // Disable writing to the depth buffer
            gl.useProgram(this.skyboxShader.program);

            // NEW: Disable attributes from the main shader that we don't use
            gl.disableVertexAttribArray(this.shader.locations.color);
            gl.disableVertexAttribArray(this.shader.locations.texcoord);
            gl.disableVertexAttribArray(this.shader.locations.normal);
            // We leave 'position' enabled as both shaders use it

            // NEW: Create a view matrix for the skybox that includes the mouse rotation
            // this.environment.modelMatrix holds the rotation from the mouse
            const skyboxViewMatrix = LIBS.multiply(this.viewMatrix, this.environment.modelMatrix);

            // MODIFIED: Use the new skyboxViewMatrix instead of this.viewMatrix
            gl.uniformMatrix4fv(this.skyboxShader.locations.Pmatrix, false, this.projMatrix);
            gl.uniformMatrix4fv(this.skyboxShader.locations.Vmatrix, false, skyboxViewMatrix);

            gl.activeTexture(gl.TEXTURE0);
            gl.bindTexture(gl.TEXTURE_CUBE_MAP, this.skyboxTexture);
            gl.uniform1i(this.skyboxShader.locations.u_skybox, 0);

            // Bind buffers
            gl.bindBuffer(gl.ARRAY_BUFFER, this.skyboxBuffers.vertex);
            // Tell WebGL how to read the position-only buffer
            gl.vertexAttribPointer(this.skyboxShader.locations.position, 3, gl.FLOAT, false, 0, 0);

            gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.skyboxBuffers.faces);

            // Draw
            gl.drawElements(gl.TRIANGLES, this.skyboxBuffers.faces_length, gl.UNSIGNED_SHORT, 0);

            gl.depthMask(true); // Re-enable depth writing

            // === 2. DRAW MAIN SCENE ===
            gl.useProgram(this.shader.program); // Switch back to the main shader

            // NEW: Re-enable attributes for the main shader
            gl.enableVertexAttribArray(this.shader.locations.color);
            gl.enableVertexAttribArray(this.shader.locations.texcoord);
            gl.enableVertexAttribArray(this.shader.locations.normal);

            // Set main scene uniforms
            gl.uniformMatrix4fv(this.shader.locations.Pmatrix, false, this.projMatrix);
            gl.uniformMatrix4fv(this.shader.locations.Vmatrix, false, this.viewMatrix);
            gl.uniform3fv(this.shader.locations.u_lightPosition, [5, 15, 10]);

            // Get the final animated matrix of the ice island
            const iceParentMatrix = this.environment.getIceIslandWorldMatrix();

            // Draw Piplup
            // First, combine the ice matrix with Piplup's static offset
            const piplupParentMatrix = LIBS.multiply(this.piplupModelMatrix, iceParentMatrix);
            // Then, draw Piplup using this final matrix
            this.piplup.draw(this.shader, piplupParentMatrix);

            const prinplupParentMatrix = LIBS.multiply(this.prinplupModelMatrix, iceParentMatrix);
            this.Prinplup.draw(this.shader, prinplupParentMatrix);

            // Draw Empoleon
            const empoleonParentMatrix = LIBS.multiply(this.empoleonModelMatrix, iceParentMatrix);
            this.empoleon.draw(this.shader, empoleonParentMatrix);

            this.environment.draw(this.shader);

            requestAnimationFrame(render);
        };
        render();
    }
}

window.addEventListener('load', () => {
    // Make sure your canvas ID matches
    new Renderer('prinplup-canvas');
});