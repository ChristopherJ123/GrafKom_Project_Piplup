import { Prinplup } from "./Resource/Prinplup.js";
import { LIBS } from "./Resource/Libs.js";
// import { Environment } from "./Environment.js";
import { ModelNode } from "./Resource/ModelNode.js";
import { Geometry } from "./Resource/Geometries.js";

// Environment khusus for my Prinplup (Michelle)
class Environment {
    constructor(gl, renderer) {
        this.gl = gl;
        this.renderer = renderer;

        this.rootNode = new ModelNode(gl); // Root of the environment scene
        this.modelMatrix = LIBS.get_I4();  // Static matrix for the whole env

        this.animatedNodes = {
            iceIsland: null
        };
        this.animationTime = 0;

        this.initParts();
    }

    initParts() {
        const gl = this.gl;

        const createTransform = (x, y, z) => {
            const m = LIBS.get_I4();
            LIBS.translateX(m, x);
            LIBS.translateY(m, y);
            LIBS.translateZ(m, z);
            return m;
        };

        const C = {
            SNOW_WHITE: [0.95, 0.98, 1.0],
            ICE_BLUE: [0.6, 0.8, 0.95],
            WATER_BLUE: [0.192, 0.502, 0.647]
        };

        // Ice Island
        const iceIslandNode = new ModelNode(gl, Geometry.generateIrregularExtrudedPolygon(6, 6, 1, C.SNOW_WHITE, 0.5));
        iceIslandNode.setBaseTransform(createTransform(0, -2.7, 0));
        this.rootNode.addChild(iceIslandNode);
        this.animatedNodes.iceIsland = iceIslandNode;

        // Water
        const waterNode = new ModelNode(gl, Geometry.generateWaterPlane(25, 25, 1, 1, C.WATER_BLUE));
        waterNode.setBaseTransform(createTransform(0, -2.6, 0));
        this.rootNode.addChild(waterNode);
    }

    // Animation update
    updateAnimation() {
        this.animationTime += 0.02;

        // Calculate float animation
        const amplitude = 0.1;
        const floatY = Math.sin(this.animationTime) * amplitude;

        // Apply animation to the ice island (combine with the base matrix)
        const T_float = LIBS.get_I4();
        LIBS.translateY(T_float, floatY);

        const iceNode = this.animatedNodes.iceIsland;
        const finalIceMatrix = LIBS.multiply(iceNode.baseMatrix, T_float);
        iceNode.setLocalTransform(finalIceMatrix);
    }

    // Ice Island: return the final computed world matrix of the island
    getIceIslandWorldMatrix() {
        this.rootNode.updateWorldMatrix(this.modelMatrix);
        return this.animatedNodes.iceIsland.worldMatrix;
    }

    draw(shader) {
        this.rootNode.updateWorldMatrix(this.modelMatrix);
        this.rootNode.draw(shader);
    }
}

class Renderer {
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        this.canvas.width = window.innerWidth;
        this.canvas.height = window.innerHeight;

        this.gl = this.canvas.getContext("webgl", { antialias: true });
        if (!this.gl) throw new Error("WebGL not supported");

        this.shader = this.createShaderProgram();
        this.Prinplup = new Prinplup(this.gl, this);
        this.environment = new Environment(this.gl, this);

        this.viewMatrix = LIBS.get_I4();
        LIBS.translateZ(this.viewMatrix, -14);
        this.projMatrix = LIBS.get_projection(30, this.canvas.width / this.canvas.height, 1, 100);
        this.animationTime = 0;

        this.audioElement = document.getElementById('backgroundMusic');

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

    loadTexture(url) {
        const gl = this.gl;
        const texture = gl.createTexture();
        gl.bindTexture(gl.TEXTURE_2D, texture);

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
        let dX = 0, dY = 0;
        let THETA = 0, PHI = 0;
        const FRICTION = 0.15;

        this.canvas.onmousedown = (e) => { drag = true; x_prev = e.pageX; y_prev = e.pageY; };
        this.canvas.onmouseup = () => { drag = false; };
        this.canvas.onmouseout = () => { drag = false; };
        this.canvas.onmousemove = (e) => {
            if (!drag) return;
            dX = (e.pageX - x_prev) * 2 * Math.PI / this.canvas.width;
            dY = (e.pageY - y_prev) * 2 * Math.PI / this.canvas.height;
            THETA += dX;
            PHI += dY;
            x_prev = e.pageX;
            y_prev = e.pageY;
        };

        this.updateRotation = () => {
            if (!drag) {
                dX *= (1 - FRICTION);
                dY *= (1 - FRICTION);
                THETA += dX;
                PHI += dY;
            }
            const rotationMatrix = LIBS.get_I4();
            LIBS.rotateY(rotationMatrix, THETA);
            LIBS.rotateX(rotationMatrix, PHI);

            // rotation applied to Prinplup's root node
            this.Prinplup.modelMatrix = rotationMatrix;
        };

        const wakeToggleButton = document.getElementById("wakeToggle");
        if (wakeToggleButton && this.audioElement) {
            wakeToggleButton.textContent = this.isAwake ? '🔇' : '🔊';

            wakeToggleButton.addEventListener('click', () => {
                this.isAwake = !this.isAwake;
                console.log("isAwake:", this.isAwake);
                if (!this.isAwake) {
                    this.audioElement.play();
                    this.audioElement.loop = true;
                    this.audioElement.currentTime = 0.5;
                    wakeToggleButton.textContent = '🔊';
                } else {
                    this.audioElement.pause();
                    this.audioElement.loop = false;
                    this.audioElement.currentTime = 0.5; // off kalau mau tidak ngulang
                    wakeToggleButton.textContent = '🔇';
                }
            });
        } else {
            if (!wakeToggleButton) console.error("Button with ID 'wakeToggle' not found!");
        }
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
            this.updateRotation();

            const timeInSeconds = now * 0.0008;
            const bodyBreathSpeed = 1.5;
            const bodyBreathAmount = 0.02;
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

            const flapSpeed = 1;
            const flapAmount = 0.2; // radians
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
            this.Prinplup.updateAnimation(animationValues);


            gl.viewport(0, 0, this.canvas.width, this.canvas.height);
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

            gl.uniformMatrix4fv(this.shader.locations.Pmatrix, false, this.projMatrix);
            gl.uniformMatrix4fv(this.shader.locations.Vmatrix, false, this.viewMatrix);

            // set the light position
            gl.uniform3fv(this.shader.locations.u_lightPosition, [5, 15, 10]);

            // 2. Get the final animated matrix of the ice island
            const iceParentMatrix = this.environment.getIceIslandWorldMatrix();

            // 3. Draw the models
            this.Prinplup.draw(this.shader, iceParentMatrix, this.isAwake);
            this.environment.draw(this.shader);

            requestAnimationFrame(render);
        };
        render();
    }
}

window.addEventListener('load', () => {
    new Renderer('prinplup-canvas');
});