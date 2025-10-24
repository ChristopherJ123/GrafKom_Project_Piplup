import { Geometry } from "./Geometries.js";
import { ModelNode } from "./ModelNode.js";
import { LIBS } from "./Libs.js";

export class Piplup {
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
        const headNode = new ModelNode(gl, Geometry.generateSphere(0.8, 0.8, 0.8, 20, 20, C.HEAD, 0.25), headTexture);
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
        LIBS.translateY(T_body, bobAmount + 1.6);

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
        this.animatedNodes.bottomBeak.setLocalTransform(bottomBeakFinal);    
    }


    draw(shader, parentMatrix) {
        // 1. Update the entire tree's matrices
        // The parentMatrix (from the ice island) is now passed in
        this.rootNode.updateWorldMatrix(parentMatrix);

        // 2. Start the recursive draw
        this.rootNode.draw(shader);
    }
}