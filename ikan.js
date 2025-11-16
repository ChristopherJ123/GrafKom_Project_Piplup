import { LIBS } from "./Resource/Libs.js";

// --- UTILITY FOR GEOMETRY ---
const Geometry = {
    generateSphere: function (a, b, c, stack, step, color) {
        const vertices = [];
        const faces = [];
        for (let i = 0; i <= stack; i++) {
            for (let j = 0; j <= step; j++) {
                const u = i / stack;
                const v = j / step;
                const theta = u * Math.PI;
                const phi = v * 2 * Math.PI;

                const x = a * Math.sin(theta) * Math.cos(phi);
                const y = b * Math.cos(theta);
                const z = c * Math.sin(theta) * Math.sin(phi);


                const texU = (v + 0.25);
                vertices.push(x, y, z, color[0], color[1], color[2], texU, u - 1.0);
            }
        }
        for (let i = 0; i < stack; i++) {
            for (let j = 0; j < step; j++) {
                const p1 = i * (step + 1) + j;
                const p2 = p1 + 1;
                const p3 = p1 + (step + 1);
                const p4 = p3 + 1;
                faces.push(p1, p2, p4, p1, p4, p3);
            }
        }
        return { vertices, faces };
    },
};

class ModelNode {
    constructor(gl, geometry = null, texture = null) {
        this.gl = gl;
        this.geometry = geometry; // The drawable geometry
        this.texture = texture;
        this.buffers = null;

        this.localMatrix = LIBS.get_I4();  // Transformation relative to its parent
        this.worldMatrix = LIBS.get_I4();  // Final transformation in world space
        this.children = [];
        this.parent = null;

        if (this.geometry) {
            this.buffers = this.createBuffers();
        }
    }

    // Function to add a child node
    addChild(node) {
        node.parent = this;
        this.children.push(node);
    }

    // Buffer creation logic
    createBuffers() {
        const vertexBuffer = this.gl.createBuffer();
        this.gl.bindBuffer(this.gl.ARRAY_BUFFER, vertexBuffer);
        this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.geometry.vertices), this.gl.STATIC_DRAW);

        const facesBuffer = this.gl.createBuffer();
        this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, facesBuffer);
        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(this.geometry.faces), this.gl.STATIC_DRAW);

        return { vertex: vertexBuffer, faces: facesBuffer, faces_length: this.geometry.faces.length };
    }

    // Set this node's local transform (e.g., the flipper's position relative to the body)
    setLocalTransform(matrix) {
        this.localMatrix = matrix;
    }

    // Recursive function to update all matrices in the tree
    updateWorldMatrix(parentWorldMatrix) {
        // Calculate our own world matrix by multiplying our local matrix
        // with our parent's world matrix.
        if (parentWorldMatrix) {
            this.worldMatrix = LIBS.multiply(this.localMatrix, parentWorldMatrix);
        } else {
            // If no parent, our world matrix is just our local matrix
            // (This is for the root node)
            this.worldMatrix = this.localMatrix;
        }

        // Now, recursively update all children
        for (const child of this.children) {
            child.updateWorldMatrix(this.worldMatrix);
        }
    }

    // Recursive function to draw this node and all its children
    draw(shader) {
        // Draw ourself (if we have geometry)
        if (this.buffers) {
            const gl = this.gl;
            gl.uniformMatrix4fv(shader.locations.Mmatrix, false, this.worldMatrix); // Use the final worldMatrix

            // (This is the same draw logic from PiplupPart)
            gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.vertex);
            gl.vertexAttribPointer(shader.locations.position, 3, gl.FLOAT, false, 4 * 8, 0);
            gl.vertexAttribPointer(shader.locations.color, 3, gl.FLOAT, false, 4 * 8, 3 * 4);
            gl.vertexAttribPointer(shader.locations.texcoord, 2, gl.FLOAT, false, 4 * 8, 6 * 4);

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

        // Now, recursively draw all children
        for (const child of this.children) {
            child.draw(shader);
        }
    }
}

// --- PIPLUP CONTAINER CLASS ---
// Manages all the parts that make up the Piplup.
class Piplup {
    constructor(gl, renderer) {
        this.gl = gl;
        this.renderer = renderer;

        // This is the single root of our scene graph
        this.rootNode = new ModelNode(gl);
        this.modelMatrix = LIBS.get_I4(); // This matrix will control the entire Piplup's rotation

        // Store references to all animated parts
        this.animatedNodes = {
            bodyNode: null,     // This node will handle the Y-translation (bob)
            bodyGeometry: null,   // This node will handle the scale (squash)
            head: null,
            topBeak: null,
            bottomBeak: null,
            leftFlipper: null,
            rightFlipper: null,
            leftLeg: null,
            rightLeg: null
        };

        // NEW: Store the base (default) transformations for animated parts
        this.baseTransforms = {
            topBeak: LIBS.get_I4(),    // NEW
            bottomBeak: LIBS.get_I4(),  // NEW
            leftFlipper: LIBS.get_I4(),
            rightFlipper: LIBS.get_I4(),
            leftLeg: LIBS.get_I4(),
            rightLeg: LIBS.get_I4(),
        };

        this.initParts();
    }

    initParts() {
        const gl = this.gl;
        // Piplup Colors
        const C = {
            BODY: [0.52, 0.80, 1.00], HEAD: [0.20, 0.38, 0.64], BEAK: [1.00, 0.84, 0.00], CAPE: [0.24, 0.42, 0.96],
            EYE_W: [1.00, 1.00, 1.00], BLACK: [0.00, 0.00, 0.00], FEET: [1.00, 0.65, 0.00], WHITE: [1.00, 1.00, 1.00]
        };

        // Helper function to create a translation matrix
        const createTransform = (x, y, z) => {
            const m = LIBS.get_I4();
            LIBS.translateX(m, x);
            LIBS.translateY(m, y);
            LIBS.translateZ(m, z);
            return m;
        };

        // Helper function
        const createOrientedTransform = (radX, radY, radZ, x, y) => {
            // Calculate z on the surface
            const termX = (x * x) / (radX * radX);
            const termY = (y * y) / (radY * radY);
            const z_on_surface = radZ * Math.sqrt(1.0 - termX - termY);

            // Add a tiny offset to prevent clipping
            const z_final = z_on_surface + 0.01;

            // Calculate rotation angles to match the surface curve
            const angleY = Math.atan2(x, z_final);
            const angleX = -Math.atan2(y, z_final);

            // Build the transformation matrix
            const m = LIBS.get_I4();
            LIBS.rotateY(m, angleY);
            LIBS.rotateX(m, angleX);

            // Apply the final position
            LIBS.set_position(m, x, y, z_final);

            return m; // Return the finished matrix
        };


        // --- Build the Hierarchy ---
        // All transforms are now RELATIVE to their parent.

        // 1. Body (Child of the root)
        // Create an empty parent node for bobbing up and down
        const bodyNode = new ModelNode(gl);
        bodyNode.setLocalTransform(createTransform(0,0,0)); // This node will be animated
        this.rootNode.addChild(bodyNode);
        this.animatedNodes.bodyNode = bodyNode; // Save for animation

        const fishTexture = this.renderer.loadTexture("Resource/fish.png");
        const headNode = new ModelNode(gl, Geometry.generateSphere(0.2, 0.8, 0.1, 20, 20, C.HEAD), fishTexture);
        const headTransform = createTransform(0, 1.5, 0); // Head is 1.5 units above the body
        headNode.setLocalTransform(headTransform);
        bodyNode.addChild(headNode);
        this.animatedNodes.head = headNode;

        // 5. Hands (Children of the Body)
        const leftHandNode = new ModelNode(gl, Geometry.generateSphere(0.2, 0.7, 0.5, 15, 15, C.BODY));
        let leftHandMatrix = createTransform(-2, 0, 0);
        // let leftHandMatrix = createTransform(0,0,0);
        LIBS.rotateZ(leftHandMatrix, LIBS.degToRad(-20));
        LIBS.rotateX(leftHandMatrix, LIBS.degToRad(-10));
        leftHandNode.setLocalTransform(leftHandMatrix);
        this.baseTransforms.leftFlipper = leftHandMatrix; // NEW: Store base pose
        bodyNode.addChild(leftHandNode);
        this.animatedNodes.leftFlipper = leftHandNode; // Save for animation
    }

    // NEW: Function to update animations
    updateAnimation(time) {
        const speed = 1; // Controls the speed of the run
        const bobAmount = Math.abs(Math.sin(time * speed) * 0.5) * 0.2 - 0.05; // Bob up and down
        const swingAngle = Math.sin(time * speed) * 0.3 - 0.4; // Swing angle in radians (~40 deg)
        const legSwingAngle = Math.sin(time * speed) * 0.5; // Legs swing a bit less

        // 1. Body Bob (with Squash and Stretch)

        // --- Calculate Squash ---
        const squashFactor = Math.abs(Math.sin(time * speed) * 0.5);
        const squashAmount = 0.1;
        const scaleY = 1.0 - squashFactor * squashAmount;
        const scaleXZ = 1.0 + squashFactor * squashAmount * 0.5;


        // 2. Left Flipper (Flaps up with top pivot)
        const T_flipper_up = LIBS.get_I4();
        // LIBS.translateY(T_flipper_up, -0.6);
        // LIBS.translateX(T_flipper_up, 0.5);

        const R_flipper = LIBS.get_I4();
        LIBS.rotateZ(R_flipper, swingAngle);

        const T_flipper_down = LIBS.get_I4();
        // LIBS.translateY(T_flipper_down, 0.6);
        // LIBS.translateX(T_flipper_down, -0.5);

        // Combine them: M_anim = T_up * (R * T_down)
        let leftFlipperAnim = LIBS.multiply(R_flipper, T_flipper_down);
        leftFlipperAnim = LIBS.multiply(T_flipper_up, leftFlipperAnim);

        // Apply animation to the base pose: M_final = M_base * M_anim
        const leftFlipperFinal = LIBS.multiply(leftFlipperAnim, this.baseTransforms.leftFlipper);
        this.animatedNodes.leftFlipper.setLocalTransform(leftFlipperFinal);
    }


    draw(shader) {
        // 1. Update the entire tree's matrices
        // this.modelMatrix is the root transform (controlled by mouse drag)
        this.rootNode.updateWorldMatrix(this.modelMatrix);

        // 2. Start the recursive draw
        this.rootNode.draw(shader);
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

        this.animationTime = 0; // NEW: Timer for animations

        this.initInputHandlers();
        this.startRenderLoop();
    }

    createShaderProgram() {
        const gl = this.gl;
        const vsSource = `
            attribute vec3 position;
            attribute vec3 color;
            attribute vec2 texcoord;
            uniform mat4 Mmatrix, Vmatrix, Pmatrix;
            varying vec3 vColor;
            varying vec2 vTexcoord;
            void main(void) {
                gl_Position = Pmatrix * Vmatrix * Mmatrix * vec4(position, 1.);
                vColor = color;
                vTexcoord = texcoord;
            }`;
        const fsSource = `
            precision mediump float;
            varying vec3 vColor;
            varying vec2 vTexcoord;
            uniform sampler2D sampler;
            uniform int u_useTexture;

            void main(void) {
                if (u_useTexture == 1) {
                    gl_FragColor = texture2D(sampler, vTexcoord);
                } else {
                    gl_FragColor = vec4(vColor, 1.);
                }
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
            Pmatrix: gl.getUniformLocation(program, "Pmatrix"),
            Vmatrix: gl.getUniformLocation(program, "Vmatrix"),
            Mmatrix: gl.getUniformLocation(program, "Mmatrix"),
            sampler: gl.getUniformLocation(program, "sampler"),
            u_useTexture: gl.getUniformLocation(program, "u_useTexture")
        };

        gl.enableVertexAttribArray(locations.position);
        gl.enableVertexAttribArray(locations.color);
        gl.enableVertexAttribArray(locations.texcoord);

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
            this.updateRotation(); // Update mouse drag rotation

            // NEW: Update the running animation
            this.animationTime += 0.05; // Increment the animation timer
            this.piplup.updateAnimation(this.animationTime);

            gl.viewport(0, 0, this.canvas.width, this.canvas.height);
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

            gl.uniformMatrix4fv(this.shader.locations.Pmatrix, false, this.projMatrix);
            gl.uniformMatrix4fv(this.shader.locations.Vmatrix, false, this.viewMatrix);

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