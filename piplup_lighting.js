import { Geometry } from "./Resource/Geometries.js";
import { LIBS } from "./Resource/Libs.js";


// --- PIPLUP PART CLASS ---
// Represents a single drawable part of the Piplup model.
class PiplupPart {
    constructor(gl, geometry, texture = null) {
        this.gl = gl;
        this.geometry = geometry;
        this.texture = texture;
        this.modelMatrix = LIBS.get_I4();
        this.buffers = this.createBuffers();
    }

    createBuffers() {
        const vertexBuffer = this.gl.createBuffer();
        this.gl.bindBuffer(this.gl.ARRAY_BUFFER, vertexBuffer);
        this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.geometry.vertices), this.gl.STATIC_DRAW);

        const facesBuffer = this.gl.createBuffer();
        this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, facesBuffer);
        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(this.geometry.faces), this.gl.STATIC_DRAW);

        return { vertex: vertexBuffer, faces: facesBuffer, faces_length: this.geometry.faces.length };
    }

    // Set the local transformation for this part (e.g., move it up, to the side, etc.)
    setTransform(transformMatrix) {
        this.modelMatrix = transformMatrix;
    }

    draw(shader, parentMatrix) {
        const gl = this.gl;
        const finalMatrix = LIBS.multiply(this.modelMatrix, parentMatrix);
        gl.uniformMatrix4fv(shader.locations.Mmatrix, false, finalMatrix);

        gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.vertex);
        const stride = 4 * (3 + 3 + 2 + 3); // 3 for pos, 3 for color, 2 for texcoord, 3 for normal
        gl.vertexAttribPointer(shader.locations.position, 3, gl.FLOAT, false, stride, 0);
        gl.vertexAttribPointer(shader.locations.color, 3, gl.FLOAT, false, stride, 3 * 4);
        gl.vertexAttribPointer(shader.locations.texcoord, 2, gl.FLOAT, false, stride, 6 * 4);
        gl.vertexAttribPointer(shader.locations.normal, 3, gl.FLOAT, false, stride, 8 * 4);


        if (this.texture) {
            gl.activeTexture(gl.TEXTURE0);
            gl.bindTexture(gl.TEXTURE_2D, this.texture);
            gl.uniform1i(shader.locations.sampler, 0);
            gl.uniform1i(shader.locations.u_useTexture, 1);
        } else {
            gl.uniform1i(shader.locations.u_useTexture, 0);
        }

        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.buffers.faces);
        gl.drawElements(gl.TRIANGLES, this.buffers.faces_length, gl.UNSIGNED_SHORT, 0);
    }
}

// --- PIPLUP CONTAINER CLASS ---
// Manages all the parts that make up the Piplup.
class Piplup {
    constructor(gl, renderer) {
        this.gl = gl;
        this.renderer = renderer;
        this.parts = [];
        this.modelMatrix = LIBS.get_I4(); // This matrix will control the entire Piplup's rotation
        this.initParts();
    }

    initParts() {
        const gl = this.gl;
        // Piplup Colors
        const C = {
            BODY: [0.52, 0.80, 1.00], HEAD: [0.20, 0.38, 0.64], BEAK: [1.00, 0.84, 0.00], CAPE: [0.24, 0.42, 0.96],
            EYE_W: [1.00, 1.00, 1.00], BLACK: [0.00, 0.00, 0.00], FEET: [1.00, 0.65, 0.00], WHITE: [1.00, 1.00, 1.00]
        };

        // Helper function to create a translation matrix using your libs.js functions
        const createTransform = (x, y, z) => {
            const m = LIBS.get_I4();
            LIBS.translateX(m, x);
            LIBS.translateY(m, y);
            LIBS.translateZ(m, z);
            return m;
        };

        const createOrientedTransform = (radX, radY, radZ, x, y) => {
            // A. Calculate the precise 'z' on the surface
            const termX = (x * x) / (radX * radX);
            const termY = (y * y) / (radY * radY);
            const z_on_surface = radZ * Math.sqrt(1.0 - termX - termY);

            // Add a tiny offset to prevent clipping
            const z_final = z_on_surface + 0.01;

            // B. Calculate rotation angles to match the surface curve
            const angleY = Math.atan2(x, z_final);
            const angleX = -Math.atan2(y, z_final);

            // C. Build the transformation matrix
            const m = LIBS.get_I4(); // Start with a fresh identity matrix

            // Apply rotations (order matters: Y-axis first, then X-axis)
            LIBS.rotateY(m, angleY);
            LIBS.rotateX(m, angleX);

            // Apply the final position
            LIBS.set_position(m, x, y, z_final);

            return m; // Return the finished matrix
        };

        const headTexture = this.renderer.loadTexture("Resource/piplup_head_texture.png");

        // Define parts and their local transformations
        const partDefinitions = [
            // Body
            { geom: Geometry.generateSphere(0.8, 0.9, 0.8, 20, 20, C.BODY), trans: LIBS.get_I4()},
            { geom: Geometry.generateCircle(0.2, 20, C.WHITE), trans: createOrientedTransform(0.8, 1.1, 0.78, -0.3, 0.25)},
            { geom: Geometry.generateCircle(0.2, 20, C.WHITE), trans: createOrientedTransform(0.8, 1.1, 0.78, 0.3, 0.25)},

            // Head
            {
                geom: Geometry.generateSphere(0.8, 0.8, 0.8, 20, 20, C.HEAD),
                trans: (() => {
                    const m = createTransform(0, 1.5, 0)
                    LIBS.rotateY(m, Math.PI / 2);
                    return m;
                })(),
                texture: headTexture
            },

            // Eyes
            { geom: Geometry.generateSphere(0.1, 0.2, 0.1, 10, 10, C.BLACK), trans: createTransform(-0.5, 1.4, 0.60)},
            { geom: Geometry.generateSphere(0.05, 0.05, 0.05, 10, 10, C.WHITE), trans: createTransform(-0.5, 1.5, 0.65)},

            { geom: Geometry.generateSphere(0.1, 0.2, 0.1, 10, 10, C.BLACK), trans: createTransform(0.5, 1.4, 0.60)},
            { geom: Geometry.generateSphere(0.05, 0.05, 0.05, 10, 10, C.WHITE), trans: createTransform(0.5, 1.5, 0.65)},

            // Beak using the new generateBeak function
            { geom: Geometry.generateBeak(0.3, 0.2, 0.6, 15, C.BEAK), trans: (() => {
                    let m = createTransform(0, 1.2, 1.2);
                    LIBS.rotateX(m, LIBS.degToRad(15));
                    return m;
                })()},
            { geom: Geometry.generateBeak(0.25, 0.15, 0.5, 15, C.BEAK), trans: (() => {
                    let m = createTransform(0, 1.15, 1.1);
                    LIBS.rotateX(m, LIBS.degToRad(5));
                    return m;
                })()},
            // Feet
            { geom: Geometry.generateSphere(0.25, 0.15, 0.5, 10, 10, C.FEET), trans: createTransform(-0.4, -1.0, 0.4)},
            { geom: Geometry.generateSphere(0.25, 0.15, 0.5, 10, 10, C.FEET), trans: createTransform(0.4, -1.0, 0.4)},
            // Legs (Body-Feet)
            { geom: Geometry.generateSphere(0.25, 0.5, 0.25, 10, 10, C.BODY), trans: createTransform(-0.4, -0.5, 0.2)},
            { geom: Geometry.generateSphere(0.25, 0.5, 0.25, 10, 10, C.BODY), trans: createTransform(0.4, -0.5, 0.2)},
            // Hands (Flippers)
            { geom: Geometry.generateSphere(0.2, 0.7, 0.5, 15, 15, C.BODY), trans: (() => {
                    let m = createTransform(-0.8, 0.1, 0.1);
                    LIBS.rotateZ(m, LIBS.degToRad(-20));
                    LIBS.rotateX(m, LIBS.degToRad(-10));
                    return m;
                })()},
            { geom: Geometry.generateSphere(0.2, 0.7, 0.5, 15, 15, C.BODY), trans: (() => {
                    let m = createTransform(0.8, 0.1, 0.1);
                    LIBS.rotateZ(m, LIBS.degToRad(20));
                    LIBS.rotateX(m, LIBS.degToRad(-10));
                    return m;
                })()},
            // Cape
            { geom: Geometry.generateTubeFromSpline(
                    // Define control points for the curve's path
                    [
                        [0.0, 0.0, -0.5],   // Start point on the lower back
                        [0.0, 0.05, -1.5],
                        [0.0, -0.2, -1.3], // First curve point
                        [0.0, -0.4, -1.8],  // Second curve point
                        [0.0, -0.7, -1.5]   // End point, slightly flared out
                    ],
                    100, // Segments for smoothness
                    0.15, // Radius of the tube (thickness of the cape)
                    20,   // Radial segments
                    [0.20, 0.38, 0.64], // Color
                ),
                // We don't need a separate transform since the points are in world space relative to the body
                trans: createTransform(0, 0, 0.5)},
            { geom: Geometry.generateSphere(0.5, 0.2, 0.4, 20, 20, C.HEAD), trans: (() => {
                    let m = createTransform(0.3, 0.8, 0.4);
                    LIBS.rotateY(m, LIBS.degToRad(-60));
                    LIBS.rotateZ(m, LIBS.degToRad(-30));
                    return m;
                })()},
            { geom: Geometry.generateSphere(0.5, 0.2, 0.4, 20, 20, C.HEAD), trans: (() => {
                    let m = createTransform(-0.3, 0.8, 0.4);
                    // LIBS.rotateX(m, LIBS.degToRad(60))
                    LIBS.rotateY(m, LIBS.degToRad(60));
                    LIBS.rotateZ(m, LIBS.degToRad(30));
                    return m;
                })()},
            // Circle cape in neck
            { geom: Geometry.generateSphere(0.7, 0.12, 0.6, 20, 20, C.HEAD), trans: (() => {
                    let m = createTransform(0, 0.8, 0);
                    // LIBS.rotateX(m, LIBS.degToRad(60))
                    // LIBS.rotateY(m, LIBS.degToRad(60));
                    // LIBS.rotateZ(m, LIBS.degToRad(30));
                    return m;
                })()},
            { geom: Geometry.generateSphere(1, 1.1, 0.25, 20, 20, C.HEAD), trans: (() => {
                    let m = createTransform(0, 0, -0.7);
                    LIBS.rotateX(m, LIBS.degToRad(20))
                    // LIBS.rotateY(m, LIBS.degToRad(60));
                    // LIBS.rotateZ(m, LIBS.degToRad(30));
                    return m;
                })()},
        ];

        partDefinitions.forEach(def => {
            const part = new PiplupPart(gl, def.geom, def.texture);
            part.setTransform(def.trans);
            this.parts.push(part);
        });
    }

    draw(shader) {
        this.parts.forEach(part => {
            part.draw(shader, this.modelMatrix);
        });
    }
}

// --- RENDERER CLASS ---
// Manages the overall WebGL scene, shaders, and render loop.
class Renderer {
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        this.canvas.width = window.innerWidth;
        this.canvas.height = window.innerHeight;

        this.gl = this.canvas.getContext("webgl", { antialias: true });
        if (!this.gl) throw new Error("WebGL not supported");

        this.shader = this.createShaderProgram();
        this.piplup = new Piplup(this.gl, this);

        this.viewMatrix = LIBS.get_I4();
        LIBS.translateZ(this.viewMatrix, -15);
        this.projMatrix = LIBS.get_projection(20, this.canvas.width / this.canvas.height, 1, 100);

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
            u_lightPosition: gl.getUniformLocation(program, "u_lightPosition")
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

        // This function is called every frame to update the Piplup's rotation
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
            this.piplup.modelMatrix = rotationMatrix;
        };
    }

    startRenderLoop() {
        const gl = this.gl;
        gl.enable(gl.DEPTH_TEST);
        gl.depthFunc(gl.LEQUAL);
        gl.clearColor(0.0, 0.0, 0.0, 0.0);
        gl.clearDepth(1.0);

        const render = () => {
            this.updateRotation();

            gl.viewport(0, 0, this.canvas.width, this.canvas.height);
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

            gl.uniformMatrix4fv(this.shader.locations.Pmatrix, false, this.projMatrix);
            gl.uniformMatrix4fv(this.shader.locations.Vmatrix, false, this.viewMatrix);

            gl.uniform3fv(this.shader.locations.u_lightPosition, [5, 15, 10]);

            this.piplup.draw(this.shader);

            requestAnimationFrame(render);
        };
        render();
    }
}

// --- START THE APPLICATION ---
window.addEventListener('load', () => {
    new Renderer('myCanvas');
});