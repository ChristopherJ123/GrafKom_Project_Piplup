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

    generateCircle: function(radius, segments, color) {
        const vertices = [];
        const faces = [];
        // Add the center vertex at (0, 0, 0)
        // Vertex format: x, y, z, r, g, b, u, v
        vertices.push(0, 0, 0, color[0], color[1], color[2], 0.5, 0.5);

        // Add vertices
        for (let i = 0; i <= segments; i++) {
            const theta = (i / segments) * 2 * Math.PI;
            const x = radius * Math.cos(theta);
            const y = radius * Math.sin(theta);

            // Texture coordinates
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

    // Elliptic paraboloid (modified radius)
    generateBeak: function(width, thickness, length, segments, color) {
        const vertices = [];
        const faces = [];

        // Tip of the beak
        vertices.push(0, 0, 0, color[0], color[1], color[2], 0, 0);

        // Build the beak with circular cross-sections
        for (let i = 1; i <= segments; i++) {
            const t = i / segments; // Parameter from 0 to 1
            const radiusScale = Math.sqrt(t); // sqrt(t) agar beak tebel di base dan tumbul di tip
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

    // Elliptic paraboloid (modified radius + cut half)
    generateBeakHalf: function(width, thickness, length, segments, color) {
        const vertices = [];
        const faces = [];

        // Tip of the beak
        vertices.push(0, 0, 0, color[0], color[1], color[2], 0, 0);

        // Build the beak with semi-circular cross-sections
        for (let i = 1; i <= segments; i++) {
            const t = i / segments;
            const radiusScale = Math.sqrt(t);
            const currentZ = -length * t;

            // Build the beak with circular cross-sections
            for (let j = 0; j <= segments; j++) {
                // Theta now goes from 0 to PI (180 deg) instead of 0 to 2*PI
                const theta = (j / segments) * Math.PI;

                const x = width * radiusScale * Math.cos(theta);
                const y = thickness * radiusScale * Math.sin(theta);
                vertices.push(x, y, currentZ, color[0], color[1], color[2], t, j / segments);
            }
        }

        // Create faces for the tip
        for (let j = 1; j <= segments; j++) {
            faces.push(0, j, j + 1);
        }

        // Create faces for the sides
        for (let i = 0; i < segments - 1; i++) {
            const ring1_start = 1 + i * (segments + 1);
            const ring2_start = 1 + (i + 1) * (segments + 1);

            for (let j = 0; j < segments; j++) {
                const p1 = ring1_start + j;
                const p2 = ring1_start + j + 1;
                const p3 = ring2_start + j;
                const p4 = ring2_start + j + 1;
                faces.push(p1, p3, p2, p2, p3, p4);
            }
        }

        // Create faces for the flat "bottom" cap (on the Y=0 plane)
        for (let i = 0; i < segments; i++) {
            const ring1_start = 1 + (i-1) * (segments + 1);
            const ring2_start = 1 + i * (segments + 1);

            // Get the four corners of the quad we are about to create
            // p1/p2 are on the previous ring, p3/p4 are on the current ring
            // The first vertex (j=0) and last (j=segments) of each ring are on the Y=0 plane.
            const p1 = (i === 0) ? 0 : ring1_start;              // First vertex of previous ring
            const p2 = (i === 0) ? 0 : ring1_start + segments;  // Last vertex of previous ring
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
                // Faktor interpolasi t (0.0-1.0)
                var t = j / segments;
                var t2 = t * t;
                var t3 = t2 * t;

                // Catmull-Rom spline formula (see google), a(tension) = 0.5
                var x = 0.5 * ((2 * p1[0]) + (-p0[0] + p2[0]) * t + (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t2 + (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t3);
                var y = 0.5 * ((2 * p1[1]) + (-p0[1] + p2[1]) * t + (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t2 + (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t3);
                var z = 0.5 * ((2 * p1[2]) + (-p0[2] + p2[2]) * t + (2 * p0[2] - 5 * p1[2] + 4 * p2[2] - p3[2]) * t2 + (-p0[2] + 3 * p1[2] - 3 * p2[2] + p3[2]) * t3);
                splinePoints.push([x, y, z]);

                // Derivative of the spline formula (the tangent(arah))
                var tx = 0.5 * ((-p0[0] + p2[0]) + 2 * (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t + 3 * (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t2);
                var ty = 0.5 * ((-p0[1] + p2[1]) + 2 * (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t + 3 * (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t2);
                var tz = 0.5 * ((-p0[2] + p2[2]) + 2 * (2 * p0[2] - 5 * p1[2] + 4 * p2[2] - p3[2]) * t + 3 * (-p0[2] + 3 * p1[2] - 3 * p2[2] + p3[2]) * t2);

                // normalisasi vektor tangent (max val 1.0)
                var mag = Math.sqrt(tx * tx + ty * ty + tz * tz);
                tangents.push([tx / mag, ty / mag, tz / mag]);
            }
        }

        // faces part
        var up = [0, 1, 0];
        for (var i = 0; i < splinePoints.length; i++) {
            var point = splinePoints[i];
            var tangent = tangents[i];

            // gimbal lock fix
            if (Math.abs(tangent[1]) > 0.999) {
                up = [1, 0, 0];
            } else {
                up = [0, 1, 0];
            }

            // Calculate a local coordinate system
            // cross product of tangent and up
            var normal = [tangent[1] * up[2] - tangent[2] * up[1], tangent[2] * up[0] - tangent[0] * up[2], tangent[0] * up[1] - tangent[1] * up[0]];
            var magN = Math.sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
            normal = [normal[0] / magN, normal[1] / magN, normal[2] / magN];

            // normalisasi normal
            // cross product of tangent and normal
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

class ModelNode {
    constructor(gl, geometry = null, texture = null) {
        this.gl = gl;
        this.geometry = geometry;
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

        // Recursively draw all children
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

        // Store the base (default) transformations for animated parts
        this.baseTransforms = {
            topBeak: LIBS.get_I4(),
            bottomBeak: LIBS.get_I4(),
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

        const headTexture = this.renderer.loadTexture("Resource/piplup_head_texture.png");

        // --- Build the Hierarchy ---
        // All transforms are now RELATIVE to their parent.

        // 1. Body (Child of the root)
        // Inviisible parent body node for bobbing up and down
        const bodyNode = new ModelNode(gl);
        bodyNode.setLocalTransform(LIBS.get_I4());
        this.rootNode.addChild(bodyNode);
        this.animatedNodes.bodyNode = bodyNode; // Save for animation

        // Body geometry as the visible part of the body (child of bodyNode)
        const bodyGeometry = new ModelNode(gl, Geometry.generateSphere(0.8, 0.9, 0.8, 20, 20, C.BODY));
        bodyGeometry.setLocalTransform(LIBS.get_I4());
        bodyNode.addChild(bodyGeometry);
        this.animatedNodes.bodyGeometry = bodyGeometry;

        // --- Body Decorations (Children of the Body) ---
        // 2a. Left White Circle
        const leftDecoCircle = new ModelNode(gl, Geometry.generateCircle(0.2, 20, C.WHITE));
        const leftDecoTransform = createOrientedTransform(0.8, 1.1, 0.78, -0.3, 0.25);
        leftDecoCircle.setLocalTransform(leftDecoTransform);
        bodyGeometry.addChild(leftDecoCircle);

        // 2b. Right White Circle
        const rightDecoCircle = new ModelNode(gl, Geometry.generateCircle(0.2, 20, C.WHITE));
        const rightDecoTransform = createOrientedTransform(0.8, 1.1, 0.78, 0.3, 0.25);
        rightDecoCircle.setLocalTransform(rightDecoTransform);
        bodyGeometry.addChild(rightDecoCircle);

        // 2. Head (Child of the Body)
        const headNode = new ModelNode(gl, Geometry.generateSphere(0.8, 0.8, 0.8, 20, 20, C.HEAD), headTexture);
        const headTransform = createTransform(0, 1.5, 0); // Head is 1.5 units above the body
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
        this.baseTransforms.leftFlipper = leftHandMatrix;
        bodyNode.addChild(leftHandNode);
        this.animatedNodes.leftFlipper = leftHandNode;

        const rightHandNode = new ModelNode(gl, Geometry.generateSphere(0.2, 0.7, 0.5, 15, 15, C.BODY));
        let rightHandMatrix = createTransform(0.8, 0.1, 0.1);
        LIBS.rotateZ(rightHandMatrix, LIBS.degToRad(20));
        LIBS.rotateX(rightHandMatrix, LIBS.degToRad(-10));
        rightHandNode.setLocalTransform(rightHandMatrix);
        this.baseTransforms.rightFlipper = rightHandMatrix;
        bodyNode.addChild(rightHandNode);
        this.animatedNodes.rightFlipper = rightHandNode;

        // 1. Left Leg
        const leftLegNode = new ModelNode(gl, Geometry.generateSphere(0.25, 0.5, 0.25, 10, 10, C.BODY));
        const leftLegTransform = createTransform(-0.4, -0.5, -0.1);
        leftLegNode.setLocalTransform(leftLegTransform);
        this.baseTransforms.leftLeg = leftLegTransform;
        bodyNode.addChild(leftLegNode);
        this.animatedNodes.leftLeg = leftLegNode;


        // 2. Left Foot (as a child of the Left Leg)
        const leftFootNode = new ModelNode(gl, Geometry.generateSphere(0.25, 0.15, 0.45, 10, 10, C.FEET));
        leftFootNode.setLocalTransform(createTransform(0, -0.5, 0.3)); // Local transform relative to leg
        leftLegNode.addChild(leftFootNode); // Add to leg

        // 3. Right Leg
        const rightLegNode = new ModelNode(gl, Geometry.generateSphere(0.25, 0.5, 0.25, 10, 10, C.BODY));
        const rightLegTransform = createTransform(0.4, -0.5, -0.1);
        rightLegNode.setLocalTransform(rightLegTransform);
        this.baseTransforms.rightLeg = rightLegTransform;
        bodyNode.addChild(rightLegNode);
        this.animatedNodes.rightLeg = rightLegNode;

        // 4. Right Foot (as a child of the Right Leg)
        const rightFootNode = new ModelNode(gl, Geometry.generateSphere(0.25, 0.15, 0.45, 10, 10, C.FEET));
        rightFootNode.setLocalTransform(createTransform(0, -0.5, 0.3)); // Local transform relative to leg
        rightLegNode.addChild(rightFootNode); // Add to leg

        // 1. The Spline Tube (Main Cape/Tail)
        const capeSplineNode = new ModelNode(gl, Geometry.generateTubeFromSpline(
            [
                [0.0, 0.0, -0.5],
                [0.0, 0.05, -1.5],
                [0.0, -0.2, -1.3],
                [0.0, -0.4, -1.8],
                [0.0, -0.7, -1.5]
            ],
            100,
            0.15,
            20,
            [0.20, 0.38, 0.64],
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
        capeNeckRing.setLocalTransform(m_neckRing);
        bodyNode.addChild(capeNeckRing);

        // 5. Back Sphere
        const capeBack = new ModelNode(gl, Geometry.generateSphere(1, 1.1, 0.25, 20, 20, C.HEAD));
        let m_back = createTransform(0, 0, -0.7);
        LIBS.rotateX(m_back, LIBS.degToRad(20));
        capeBack.setLocalTransform(m_back);
        bodyNode.addChild(capeBack);
    }

    // Function to update animations
    updateAnimation(time) {
        const speed = 2; // Controls the speed of the run
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
        // Apply Translation to the parent bobber
        this.animatedNodes.bodyNode.setLocalTransform(T_body);

        // Apply Scale ONLY to the body geometry
        this.animatedNodes.bodyGeometry.setLocalTransform(S_body);


        // 2. Left Flipper (Flaps up with top pivot)
        const T_flipper_up = LIBS.get_I4();
        LIBS.translateX(T_flipper_up, 0.5);
        LIBS.translateY(T_flipper_up, -0.6);

        let R_flipper = LIBS.get_I4();
        let swingAxis = [0, 0, 1.0];
        LIBS.rotateAroundAxis(R_flipper, swingAxis, swingAngle);

        const T_flipper_down = LIBS.get_I4();
        LIBS.translateX(T_flipper_down, -0.5);
        LIBS.translateY(T_flipper_down, 0.6);


        // M_anim = T_up * (R * T_down)
        let leftFlipperAnim = LIBS.multiply(R_flipper, T_flipper_down);
        leftFlipperAnim = LIBS.multiply(T_flipper_up, leftFlipperAnim);

        // Base: M_final = M_base * M_anim
        const leftFlipperFinal = LIBS.multiply(this.baseTransforms.leftFlipper, leftFlipperAnim);
        this.animatedNodes.leftFlipper.setLocalTransform(leftFlipperFinal);


        // 3. Right Flipper
        const T_flipper_up_right = LIBS.get_I4();
        LIBS.translateY(T_flipper_up_right, -0.6);
        LIBS.translateX(T_flipper_up_right, -0.5); // Mirrored X

        const R_flipper_right = LIBS.get_I4();
        swingAxis = [0, 0, 1.0];
        LIBS.rotateAroundAxis(R_flipper_right, swingAxis, swingAngle + 0.8);

        const T_flipper_down_right = LIBS.get_I4();
        LIBS.translateY(T_flipper_down_right, 0.6);
        LIBS.translateX(T_flipper_down_right, 0.5); // Mirrored X

        // Use the right-flipper matrices
        // M_anim = T_up_right * (R_right * T_down_right)
        let rightFlipperAnim = LIBS.multiply(R_flipper_right, T_flipper_down_right);
        rightFlipperAnim = LIBS.multiply(T_flipper_up_right, rightFlipperAnim);

        // Base: M_final = M_base * M_anim
        const rightFlipperFinal = LIBS.multiply(this.baseTransforms.rightFlipper, rightFlipperAnim);
        this.animatedNodes.rightFlipper.setLocalTransform(rightFlipperFinal);

        // 4. Left Leg (Swings backward)
        const leftLegAnim = LIBS.get_I4();
        LIBS.rotateX(leftLegAnim, -legSwingAngle);
        const leftLegFinal = LIBS.multiply(leftLegAnim, this.baseTransforms.leftLeg);
        this.animatedNodes.leftLeg.setLocalTransform(leftLegFinal);

        // 5. Right Leg (Swings forward)
        const rightLegAnim = LIBS.get_I4();
        LIBS.rotateX(rightLegAnim, legSwingAngle);
        const rightLegFinal = LIBS.multiply(rightLegAnim, this.baseTransforms.rightLeg);
        this.animatedNodes.rightLeg.setLocalTransform(rightLegFinal);


        // 6. Beak Animation
        const beakAngle = (Math.sin(time * speed * 1.5 + 2) + 1) * 0.5 * 0.2;

        // Pivot hinge di z = -0.3 relatif dari beak
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