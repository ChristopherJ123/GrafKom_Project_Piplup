import { LIBS } from "./Resource/Libs.js";
// --- UTILITY FOR GEOMETRY ---
// We place the geometry generation logic in its own object to keep things organized.
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


                // We add 0.25 to shift the coordinates, so v=0.25 (front) maps to U=0.5 (texture center)
                const texU = (v + 0.25); // Use modulo to wrap around
                vertices.push(x, y, z, color[0], color[1], color[2], texU, u - 1.0);            }
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

    /**
     * Generates a blunted, beak-like shape.
     * @param {number} width - Max width (x-axis).
     * @param {number} thickness - Max thickness (y-axis).
     * @param {number} length - The length of the beak (z-axis).
     * @param {number} segments - The number of segments for resolution.
     * @param {Array<number>} color - The RGB color array.
     * @returns {{vertices: Array<number>, faces: Array<number>}}
     */
    generateCircle: function(radius, segments, color) {
        const vertices = [];
        const faces = [];
        // Add the center vertex at (0, 0, 0)
        // Vertex format: x, y, z, r, g, b, u, v
        vertices.push(0, 0, 0, color[0], color[1], color[2], 0.5, 0.5);

        // Add vertices for the circumference
        for (let i = 0; i <= segments; i++) {
            const theta = (i / segments) * 2 * Math.PI;
            const x = radius * Math.cos(theta);
            const y = radius * Math.sin(theta);

            // Texture coordinates (map circular coords to square UV space)
            const u = (x / radius + 1) / 2;
            const v = (y / radius + 1) / 2;

            vertices.push(x, y, 0, color[0], color[1], color[2], u, v);
        }
        // Create the faces using the center vertex (index 0)
        for (let i = 1; i <= segments; i++) {
            const centerIndex = 0;
            const currentIndex = i;
            const nextIndex = i + 1;
            faces.push(centerIndex, currentIndex, nextIndex);
        }
        return { vertices, faces };
    },

    // Elliptic paraboloid
    generateBeak: function(width, thickness, length, segments, color) {
        const vertices = [];
        const faces = [];

        // Tip of the beak
        vertices.push(0, 0, 0, color[0], color[1], color[2], 0, 0);

        // Build the beak with circular cross-sections
        for (let i = 1; i <= segments; i++) {

            // Key functions
            const t = i / segments; // Parameter from 0 to 1
            const radiusScale = Math.sqrt(t); // Use sqrt(t) to make the beak fatter at the base and blunter at the tip
            const currentZ = -length * t; // Move along the negative Z-axis

            for (let j = 0; j < segments; j++) {
                const theta = (j / segments) * 2 * Math.PI;
                const x = width * radiusScale * Math.cos(theta);
                const y = thickness * radiusScale * Math.sin(theta);
                vertices.push(x, y, currentZ, color[0], color[1], color[2], t, j / segments);
            }
        }

        // Create faces for the tip
        for (let j = 1; j <= segments; j++) {
            faces.push(0, j, (j % segments) + 1);
        }

        // Create faces for the sides
        for (let i = 0; i < segments - 1; i++) {
            const ring1_start = 1 + i * segments;
            const ring2_start = 1 + (i + 1) * segments;
            for (let j = 0; j < segments; j++) {
                const p1 = ring1_start + j;
                const p2 = ring1_start + ((j + 1) % segments);
                const p3 = ring2_start + j;
                const p4 = ring2_start + ((j + 1) % segments);
                faces.push(p1, p3, p2, p2, p3, p4);
            }
        }

        return { vertices, faces };
    },


    /**
     * Generates one half (the top half, in positive Y) of a blunted, beak-like shape.
     * @param {number} width - Max width (x-axis).
     * @param {number} thickness - Max thickness (y-axis).
     * @param {number} length - The length of the beak (z-axis).
     * @param {number} segments - The number of segments for resolution.
     * @param {Array<number>} color - The RGB color array.
     * @returns {{vertices: Array<number>, faces: Array<number>}}
     */
    generateBeakHalf: function(width, thickness, length, segments, color) {
        const vertices = [];
        const faces = [];

        // Tip of the beak (remains the same)
        vertices.push(0, 0, 0, color[0], color[1], color[2], 0, 0);

        // Build the beak with semi-circular cross-sections
        for (let i = 1; i <= segments; i++) {
            // Key functions (remain the same)
            const t = i / segments;
            const radiusScale = Math.sqrt(t);
            const currentZ = -length * t;

            // MODIFIED: Loop from 0 to segments (inclusive) to create segments+1 vertices
            for (let j = 0; j <= segments; j++) {
                // MODIFIED: Theta now goes from 0 to PI (180 deg) instead of 0 to 2*PI
                const theta = (j / segments) * Math.PI;

                const x = width * radiusScale * Math.cos(theta);
                const y = thickness * radiusScale * Math.sin(theta); // This will only be >= 0
                vertices.push(x, y, currentZ, color[0], color[1], color[2], t, j / segments);
            }
        }

        // Create faces for the curved tip
        // MODIFIED: Loop 'j' from 1 to segments (not <=), no modulo needed
        for (let j = 1; j <= segments; j++) {
            // We connect (0) -> (j) -> (j+1)
            faces.push(0, j, j + 1);
        }

        // Create faces for the curved sides
        for (let i = 0; i < segments - 1; i++) {
            // MODIFIED: We have (segments + 1) vertices per ring
            const ring1_start = 1 + i * (segments + 1);
            const ring2_start = 1 + (i + 1) * (segments + 1);

            // MODIFIED: Loop 'j' from 0 to segments-1 (to create 'segments' quads)
            for (let j = 0; j < segments; j++) {
                // MODIFIED: No modulo needed, just use j and j+1
                const p1 = ring1_start + j;
                const p2 = ring1_start + j + 1;
                const p3 = ring2_start + j;
                const p4 = ring2_start + j + 1;
                faces.push(p1, p3, p2, p2, p3, p4);
            }
        }

        // NEW: Create faces for the flat "bottom" cap (on the Y=0 plane)
        // This connects all the edge vertices to seal the model
        for (let i = 0; i < segments; i++) {
            const ring1_start = 1 + (i-1) * (segments + 1);
            const ring2_start = 1 + i * (segments + 1);

            // Get the four corners of the quad we are about to create
            // p1/p2 are on the previous ring, p3/p4 are on the current ring
            // The first vertex (j=0) and last (j=segments) of each ring are on the Y=0 plane.
            const p1 = (i === 0) ? 0 : ring1_start;              // First vertex of previous ring (or tip)
            const p2 = (i === 0) ? 0 : ring1_start + segments;  // Last vertex of previous ring (or tip)
            const p3 = ring2_start;                             // First vertex of current ring
            const p4 = ring2_start + segments;                  // Last vertex of current ring

            if (i === 0) {
                // Create the first triangle at the tip
                faces.push(p1, p4, p3);
            } else {
                // Create the quads for the rest of the cap
                faces.push(p1, p2, p3, p2, p4, p3);
            }
        }

        return { vertices, faces };
    },

    generateTubeFromSpline: function(controlPoints, segments, radius, radialSegments, color) {
        var vertices = [];
        var faces = [];
        var splinePoints = [];
        var tangents = [];

        var points = [];
        points.push(controlPoints[0]);
        controlPoints.forEach(p => points.push(p));
        points.push(controlPoints[controlPoints.length - 1]);

        for (var i = 1; i < points.length - 2; i++) {
            var p0 = points[i - 1];
            var p1 = points[i];
            var p2 = points[i + 1];
            var p3 = points[i + 2];

            for (var j = 0; j <= segments; j++) {
                var t = j / segments;
                var t2 = t * t;
                var t3 = t2 * t;

                var x = 0.5 * ((2 * p1[0]) + (-p0[0] + p2[0]) * t + (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t2 + (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t3);
                var y = 0.5 * ((2 * p1[1]) + (-p0[1] + p2[1]) * t + (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t2 + (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t3);
                var z = 0.5 * ((2 * p1[2]) + (-p0[2] + p2[2]) * t + (2 * p0[2] - 5 * p1[2] + 4 * p2[2] - p3[2]) * t2 + (-p0[2] + 3 * p1[2] - 3 * p2[2] + p3[2]) * t3);
                splinePoints.push([x, y, z]);

                var tx = 0.5 * ((-p0[0] + p2[0]) + 2 * (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t + 3 * (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t2);
                var ty = 0.5 * ((-p0[1] + p2[1]) + 2 * (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t + 3 * (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t2);
                var tz = 0.5 * ((-p0[2] + p2[2]) + 2 * (2 * p0[2] - 5 * p1[2] + 4 * p2[2] - p3[2]) * t + 3 * (-p0[2] + 3 * p1[2] - 3 * p2[2] + p3[2]) * t2);

                var mag = Math.sqrt(tx * tx + ty * ty + tz * tz);
                tangents.push([tx / mag, ty / mag, tz / mag]);
            }
        }

        var up = [0, 1, 0];
        for (var i = 0; i < splinePoints.length; i++) {
            var point = splinePoints[i];
            var tangent = tangents[i];

            if (Math.abs(tangent[1]) > 0.999) {
                up = [1, 0, 0];
            } else {
                up = [0, 1, 0];
            }

            var normal = [tangent[1] * up[2] - tangent[2] * up[1], tangent[2] * up[0] - tangent[0] * up[2], tangent[0] * up[1] - tangent[1] * up[0]];
            var magN = Math.sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
            normal = [normal[0] / magN, normal[1] / magN, normal[2] / magN];

            var binormal = [tangent[1] * normal[2] - tangent[2] * normal[1], tangent[2] * normal[0] - tangent[0] * normal[2], tangent[0] * normal[1] - tangent[1] * normal[0]];

            for (var j = 0; j <= radialSegments; j++) {
                var theta = (j / radialSegments) * 2 * Math.PI;
                var x = point[0] + radius * (Math.cos(theta) * normal[0] + Math.sin(theta) * binormal[0]);
                var y = point[1] + radius * (Math.cos(theta) * normal[1] + Math.sin(theta) * binormal[1]);
                var z = point[2] + radius * (Math.cos(theta) * normal[2] + Math.sin(theta) * binormal[2]);
                vertices.push(x, y, z);

                // Color
                vertices.push(color[0], color[1], color[2], j / radialSegments, i / splinePoints.length);
            }
        }

        for (var i = 0; i < splinePoints.length - 1; i++) {
            for (var j = 0; j < radialSegments; j++) {
                var p1 = i * (radialSegments + 1) + j;
                var p2 = p1 + 1;
                var p3 = (i + 1) * (radialSegments + 1) + j;
                var p4 = p3 + 1;
                faces.push(p1, p2, p4);
                faces.push(p1, p4, p3);
            }
        }
        return { vertices, faces };
    }
};

// This class will replace BOTH Piplup and PiplupPart
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

    // (This is the same buffer creation logic from PiplupPart)
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

    // NEW: Recursive function to update all matrices in the tree
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

    // NEW: Recursive function to draw this node and all its children
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

        // MODIFIED: Store references to all animated parts
        this.animatedNodes = {
            bodyNode: null,     // NEW: This node will handle the Y-translation (bob)
            bodyGeometry: null,   // NEW: This node will handle the scale (squash)
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

        // --- Build the Hierarchy ---
        // All transforms are now RELATIVE to their parent.

// 1. Body (Child of the root)
        // MODIFIED: Create an empty parent node for bobbing
        const bodyNode = new ModelNode(gl);
        bodyNode.setLocalTransform(LIBS.get_I4()); // This node will be animated
        this.rootNode.addChild(bodyNode);
        this.animatedNodes.bodyNode = bodyNode; // Save for animation

        // NEW: Create the visible body geometry as a child of the bobber
        const bodyGeometry = new ModelNode(gl, Geometry.generateSphere(0.8, 0.9, 0.8, 20, 20, C.BODY));
        bodyGeometry.setLocalTransform(LIBS.get_I4()); // This node will be scaled
        bodyNode.addChild(bodyGeometry); // Add to the bobber
        this.animatedNodes.bodyGeometry = bodyGeometry; // Save for animation

        // --- Body Decorations (Children of the Body) ---
        // (Assuming 'bodyNode' is your main body ModelNode)

        // 1. Left White Circle
        const leftDecoCircle = new ModelNode(gl, Geometry.generateCircle(0.2, 20, C.WHITE));
        // We use the original function to get the correct transformation matrix
        const leftDecoTransform = createOrientedTransform(0.8, 1.1, 0.78, -0.3, 0.25);
        leftDecoCircle.setLocalTransform(leftDecoTransform);
        bodyGeometry.addChild(leftDecoCircle);

        // 2. Right White Circle
        const rightDecoCircle = new ModelNode(gl, Geometry.generateCircle(0.2, 20, C.WHITE));
        const rightDecoTransform = createOrientedTransform(0.8, 1.1, 0.78, 0.3, 0.25);
        rightDecoCircle.setLocalTransform(rightDecoTransform);
        bodyGeometry.addChild(rightDecoCircle);

        // 2. Head (Child of the Body)
        const headNode = new ModelNode(gl, Geometry.generateSphere(0.8, 0.8, 0.8, 20, 20, C.HEAD), headTexture);
        const headTransform = createTransform(0, 1.5, 0); // Head is 1.5 units *above the body*
        headNode.setLocalTransform(headTransform);
        bodyNode.addChild(headNode);
        this.animatedNodes.head = headNode;

        // 3. Eyes (Children of the Head)
        // Eye position is now relative to the Head (at 0, 1.5, 0)
        // Original eye: (-0.5, 1.4, 0.60)
        // New relative eye: (-0.5, 1.4 - 1.5, 0.60) -> (-0.5, -0.1, 0.60)
        const leftEyeNode = new ModelNode(gl, Geometry.generateSphere(0.1, 0.2, 0.1, 10, 10, C.BLACK));
        leftEyeNode.setLocalTransform(createTransform(-0.5, -0.1, 0.60));
        headNode.addChild(leftEyeNode);

        const leftEyeLightNode = new ModelNode(gl, Geometry.generateSphere(0.05, 0.05, 0.05, 10, 10, C.WHITE));
        leftEyeLightNode.setLocalTransform(createTransform(-0.5, 0, 0.65));
        headNode.addChild(leftEyeLightNode);

        // Right eye
        const rightEyeNode = new ModelNode(gl, Geometry.generateSphere(0.1, 0.2, 0.1, 10, 10, C.BLACK));
        rightEyeNode.setLocalTransform(createTransform(0.5, -0.1, 0.60));
        headNode.addChild(rightEyeNode);

        const rightEyeLightNode = new ModelNode(gl, Geometry.generateSphere(0.05, 0.05, 0.05, 10, 10, C.WHITE));
        rightEyeLightNode.setLocalTransform(createTransform(0.5, 0, 0.65));
        headNode.addChild(rightEyeLightNode);


        // 4. Beak (Child of the Head)
        // Original beak: (0, 1.2, 1.2)
        // New relative beak: (0, 1.2 - 1.5, 1.2) -> (0, -0.3, 1.2)
        const beakUpNode = new ModelNode(gl, Geometry.generateBeakHalf(0.25, 0.3, 0.3, 15, C.BEAK));
        let beakUpMatrix = createTransform(0, -0.35, 1);
        LIBS.rotateX(beakUpMatrix, LIBS.degToRad(15));
        beakUpNode.setLocalTransform(beakUpMatrix);
        headNode.addChild(beakUpNode);
        this.animatedNodes.topBeak = beakUpNode;
        this.baseTransforms.topBeak = beakUpMatrix;

        const beakDownNode = new ModelNode(gl, Geometry.generateBeakHalf(0.25, 0.25, 0.3, 15, C.BEAK));
        let beakDownMatrix = createTransform(0, -0.35, 1);
        LIBS.rotateX(beakDownMatrix, LIBS.degToRad(-25));
        LIBS.rotateZ(beakDownMatrix, LIBS.degToRad(180));
        beakDownNode.setLocalTransform(beakDownMatrix);
        headNode.addChild(beakDownNode);
        this.animatedNodes.bottomBeak = beakDownNode;
        this.baseTransforms.bottomBeak = beakDownMatrix;

        // 5. Hands (Children of the Body)
        const leftHandNode = new ModelNode(gl, Geometry.generateSphere(0.2, 0.7, 0.5, 15, 15, C.BODY));
        let leftHandMatrix = createTransform(-0.8, 0.1, 0.1);
        LIBS.rotateZ(leftHandMatrix, LIBS.degToRad(-20));
        LIBS.rotateX(leftHandMatrix, LIBS.degToRad(-10));
        leftHandNode.setLocalTransform(leftHandMatrix);
        this.baseTransforms.leftFlipper = leftHandMatrix; // NEW: Store base pose
        bodyNode.addChild(leftHandNode);
        this.animatedNodes.leftFlipper = leftHandNode; // Save for animation

        const rightHandNode = new ModelNode(gl, Geometry.generateSphere(0.2, 0.7, 0.5, 15, 15, C.BODY));
        let rightHandMatrix = createTransform(0.8, 0.1, 0.1);
        LIBS.rotateZ(rightHandMatrix, LIBS.degToRad(20));
        LIBS.rotateX(rightHandMatrix, LIBS.degToRad(-10));
        rightHandNode.setLocalTransform(rightHandMatrix);
        this.baseTransforms.rightFlipper = rightHandMatrix; // NEW: Store base pose
        bodyNode.addChild(rightHandNode);
        this.animatedNodes.rightFlipper = rightHandNode; // MODIFIED: Corrected typo (was leftFlipper)

        // ... (add right flipper as a child of bodyNode) ...

        // --- Cape and Back Details (Children of the Body) ---
        // (Assuming 'bodyNode' is the ModelNode for the Piplup's body)

        // --- Legs and Feet (Children of Body) ---

        // (Assuming 'bodyNode' is your main body ModelNode)
        // (And C.FEET and C.BODY are defined in your colors)

        // 1. Left Leg
        const leftLegNode = new ModelNode(gl, Geometry.generateSphere(0.25, 0.5, 0.25, 10, 10, C.BODY));
        const leftLegTransform = createTransform(-0.4, -0.5, -0.1); // NEW: Capture transform
        leftLegNode.setLocalTransform(leftLegTransform);
        this.baseTransforms.leftLeg = leftLegTransform; // NEW: Store base pose
        bodyNode.addChild(leftLegNode);
        this.animatedNodes.leftLeg = leftLegNode; // NEW: Save for animation


        // 2. Left Foot (as a child of the Left Leg)
        const leftFootNode = new ModelNode(gl, Geometry.generateSphere(0.25, 0.15, 0.45, 10, 10, C.FEET));
        // The transform is (0, -0.5, 0.2) *relative* to the leg's position
        leftFootNode.setLocalTransform(createTransform(0, -0.5, 0.3));
        leftLegNode.addChild(leftFootNode); // <-- Add to leg, NOT body

        // 3. Right Leg
        const rightLegNode = new ModelNode(gl, Geometry.generateSphere(0.25, 0.5, 0.25, 10, 10, C.BODY));
        const rightLegTransform = createTransform(0.4, -0.5, -0.1); // NEW: Capture transform
        rightLegNode.setLocalTransform(rightLegTransform);
        this.baseTransforms.rightLeg = rightLegTransform; // NEW: Store base pose
        bodyNode.addChild(rightLegNode);
        this.animatedNodes.rightLeg = rightLegNode; // NEW: Save for animation

        // 4. Right Foot (as a child of the Right Leg)
        const rightFootNode = new ModelNode(gl, Geometry.generateSphere(0.25, 0.15, 0.45, 10, 10, C.FEET));
        // The transform is (0, -0.5, 0.2) *relative* to the leg's position
        rightFootNode.setLocalTransform(createTransform(0, -0.5, 0.3));
        rightLegNode.addChild(rightFootNode); // <-- Add to leg, NOT body

        // 1. The Spline Tube (Main Cape/Tail)
        const capeSplineNode = new ModelNode(gl, Geometry.generateTubeFromSpline(
            [
                [0.0, 0.0, -0.5],   // Start point on the lower back
                [0.0, 0.05, -1.5],
                [0.0, -0.2, -1.3], // First curve point
                [0.0, -0.4, -1.8],  // Second curve point
                [0.0, -0.7, -1.5]   // End point, slightly flared out
            ],
            100, // Segments
            0.15, // Radius
            20,   // Radial segments
            [0.20, 0.38, 0.64] // Color
        ));
        capeSplineNode.setLocalTransform(createTransform(0, 0, 0.5));
        bodyNode.addChild(capeSplineNode);

        // 2. Right Shoulder Sphere
        const capeRightShoulder = new ModelNode(gl, Geometry.generateSphere(0.5, 0.2, 0.4, 20, 20, C.HEAD));
        let m_rightShoulder = createTransform(0.3, 0.8, 0.4);
        LIBS.rotateY(m_rightShoulder, LIBS.degToRad(-60));
        LIBS.rotateZ(m_rightShoulder, LIBS.degToRad(-30));
        capeRightShoulder.setLocalTransform(m_rightShoulder);
        bodyNode.addChild(capeRightShoulder);

        // 3. Left Shoulder Sphere
        const capeLeftShoulder = new ModelNode(gl, Geometry.generateSphere(0.5, 0.2, 0.4, 20, 20, C.HEAD));
        let m_leftShoulder = createTransform(-0.3, 0.8, 0.4);
        LIBS.rotateY(m_leftShoulder, LIBS.degToRad(60));
        LIBS.rotateZ(m_leftShoulder, LIBS.degToRad(30));
        capeLeftShoulder.setLocalTransform(m_leftShoulder);
        bodyNode.addChild(capeLeftShoulder);

        // 4. Neck Ring Sphere
        const capeNeckRing = new ModelNode(gl, Geometry.generateSphere(0.7, 0.12, 0.6, 20, 20, C.HEAD));
        let m_neckRing = createTransform(0, 0.8, 0);
        // (Original rotations were commented out, so I've left them out)
        capeNeckRing.setLocalTransform(m_neckRing);
        bodyNode.addChild(capeNeckRing);

        // 5. Back Sphere
        const capeBack = new ModelNode(gl, Geometry.generateSphere(1, 1.1, 0.25, 20, 20, C.HEAD));
        let m_back = createTransform(0, 0, -0.7);
        LIBS.rotateX(m_back, LIBS.degToRad(20));
        capeBack.setLocalTransform(m_back);
        bodyNode.addChild(capeBack);

        // ... (Continue for all other parts, attaching them to either the body or head) ...

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

        // --- Create Translation Matrix ---
        const T_body = LIBS.get_I4();
        LIBS.translateY(T_body, bobAmount);

        // --- Create Scale Matrix ---
        const S_body = LIBS.get_I4();
        LIBS.scaleY(S_body, scaleY);
        LIBS.scaleX(S_body, scaleXZ);
        LIBS.scaleZ(S_body, scaleXZ);

        // --- Apply Matrices ---
        // MODIFIED: Apply Translation to the parent bobber
        this.animatedNodes.bodyNode.setLocalTransform(T_body);

        // MODIFIED: Apply Scale ONLY to the body geometry
        this.animatedNodes.bodyGeometry.setLocalTransform(S_body);


        // 2. Left Flipper (Flaps up with top pivot)
        // Create the animation matrix: M_anim = T(0.7) * R(angle) * T(-0.7)
        const T_flipper_up = LIBS.get_I4();
        LIBS.translateY(T_flipper_up, -0.6);
        LIBS.translateX(T_flipper_up, 0.5);

        const R_flipper = LIBS.get_I4();
        // MODIFIED: Changed from rotateX to rotateZ
        LIBS.rotateZ(R_flipper, swingAngle);

        const T_flipper_down = LIBS.get_I4();
        LIBS.translateY(T_flipper_down, 0.6);
        LIBS.translateX(T_flipper_down, -0.5);

        // Combine them: M_anim = T_up * (R * T_down)
        let leftFlipperAnim = LIBS.multiply(R_flipper, T_flipper_down);
        leftFlipperAnim = LIBS.multiply(T_flipper_up, leftFlipperAnim);

        // Apply animation to the base pose: M_final = M_base * M_anim
        const leftFlipperFinal = LIBS.multiply(this.baseTransforms.leftFlipper, leftFlipperAnim);
        this.animatedNodes.leftFlipper.setLocalTransform(leftFlipperFinal);


        // 3. Right Flipper (Flaps down with top pivot)
// --- NEW: Create mirrored pivot matrices for the right flipper ---
// Pivot point is at (x: 0.5, y: 0.6) relative to its center
        const T_flipper_up_right = LIBS.get_I4();
        LIBS.translateY(T_flipper_up_right, -0.6);
        LIBS.translateX(T_flipper_up_right, -0.5); // Mirrored X

        const R_flipper_right = LIBS.get_I4();
        LIBS.rotateZ(R_flipper_right, swingAngle + 0.8); // Note the negative angle

        const T_flipper_down_right = LIBS.get_I4();
        LIBS.translateY(T_flipper_down_right, 0.6);
        LIBS.translateX(T_flipper_down_right, 0.5); // Mirrored X

// --- MODIFIED: Use the new right-flipper matrices ---
// Combine them: M_anim = T_up_right * (R_right * T_down_right)
        let rightFlipperAnim = LIBS.multiply(R_flipper_right, T_flipper_down_right);
        rightFlipperAnim = LIBS.multiply(T_flipper_up_right, rightFlipperAnim);

// Apply animation to the base pose: M_final = M_base * M_anim
        const rightFlipperFinal = LIBS.multiply(this.baseTransforms.rightFlipper, rightFlipperAnim);
        this.animatedNodes.rightFlipper.setLocalTransform(rightFlipperFinal);

        // 4. Left Leg (Swings backward, opposite of left flipper)
        const leftLegAnim = LIBS.get_I4();
        LIBS.rotateX(leftLegAnim, -legSwingAngle);
        const leftLegFinal = LIBS.multiply(leftLegAnim, this.baseTransforms.leftLeg);
        this.animatedNodes.leftLeg.setLocalTransform(leftLegFinal);

        // 5. Right Leg (Swings forward, opposite of right flipper)
        const rightLegAnim = LIBS.get_I4();
        LIBS.rotateX(rightLegAnim, legSwingAngle);
        const rightLegFinal = LIBS.multiply(rightLegAnim, this.baseTransforms.rightLeg);
        this.animatedNodes.rightLeg.setLocalTransform(rightLegFinal);


// 6. Beak "Talking" Animation
        // Creates a small "chirp" animation.
        // (Math.sin + 1)/2 gives a 0-to-1 range. * 0.2 gives a max 0.2 rad angle.
        const beakAngle = (Math.sin(time * speed * 1.5 + 2) + 1) * 0.5 * 0.2;

        // Define the pivot matrices. Hinge is at z = -0.3
        const T_beak_hinge = LIBS.get_I4();
        LIBS.translateZ(T_beak_hinge, 0.7); // 1. Move hinge to origin

        const T_beak_back = LIBS.get_I4();
        LIBS.translateZ(T_beak_back, -0.7); // 3. Move hinge back

        // --- Top Beak ---
        const R_top_beak = LIBS.get_I4();
        LIBS.rotateX(R_top_beak, -beakAngle); // 2. Rotate (opens up)

        // Combine: M_anim = T_back * R_anim * T_hinge
        let topBeakAnim = LIBS.multiply(R_top_beak, T_beak_hinge);
        topBeakAnim = LIBS.multiply(T_beak_back, topBeakAnim);

        // Apply to base pose: M_final = M_base * M_anim
        const topBeakFinal = LIBS.multiply(this.baseTransforms.topBeak, topBeakAnim);
        this.animatedNodes.topBeak.setLocalTransform(topBeakFinal);

        // --- Bottom Beak ---
        const R_bottom_beak = LIBS.get_I4();
        LIBS.rotateX(R_bottom_beak, beakAngle); // 2. Rotate (opens down)

        // Combine: M_anim = T_back * R_anim * T_hinge
        let bottomBeakAnim = LIBS.multiply(R_bottom_beak, T_beak_hinge);
        bottomBeakAnim = LIBS.multiply(T_beak_back, bottomBeakAnim);

        // Apply to base pose: M_final = M_base * M_anim
        const bottomBeakFinal = LIBS.multiply(this.baseTransforms.bottomBeak, bottomBeakAnim);
        this.animatedNodes.bottomBeak.setLocalTransform(bottomBeakFinal);    }


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