import { Geometry } from "./Geometries.js";
import { LIBS } from "./Libs.js";
import { ModelNode } from "./ModelNode.js";

export class Prinplup {
    constructor(gl, renderer) {
        this.gl = gl;
        this.renderer = renderer;

        this.rootNode = new ModelNode(gl);

        // MODIFIED: Set the static Y-offset here
        this.modelMatrix = LIBS.get_I4();
        // LIBS.translateY(this.modelMatrix, 2.5); // Prinplup's base offset

        // Store references to nodes that need animation
        this.animatedNodes = {
            breathingNode: null,
            leftLegGroup: null, 
            rightLegGroup: null,
            eyeNodes: [],
            leftHand: null,
            rightHand: null,
            breathEffects: []
        };

        // NEW: Add a place to store base poses
        this.baseTransforms = {
            leftHand: LIBS.get_I4(),
            rightHand: LIBS.get_I4(),
            leftLegGroup: LIBS.get_I4(),
            rightLegGroup: LIBS.get_I4()
        };

        this.initParts();
    }

    initParts() {
        const gl = this.gl;
        const C = {
            BODY: [0.61, 0.84, 0.89], HEAD: [0.20, 0.38, 0.64], BEAK: [1.00, 0.84, 0.00], CAPE: [0.24, 0.42, 0.96],
            EYE_W: [1.00, 1.00, 1.00], BLACK: [0.00, 0.00, 0.00], FEET: [1.00, 0.65, 0.00], WHITE: [1.00, 1.00, 1.00], TAIL: [0.08, 0.32, 0.60]
        };

        const createTransform = (x, y, z) => {
            const m = LIBS.get_I4();
            LIBS.translateX(m, x);
            LIBS.translateY(m, y);
            LIBS.translateZ(m, z);
            return m;
        };

        const createOrientedTransform = (radX, radY, radZ, x, y) => {
            const termX = (x * x) / (radX * radX);
            const termY = (y * y) / (radY * radZ); // Bug? Should be radY * radY
            const z_on_surface = radZ * Math.sqrt(1.0 - termX - termY);
            const z_final = z_on_surface + 0.01;
            const angleY = Math.atan2(x, z_final);
            const angleX = -Math.atan2(y, z_final);
            const m = LIBS.get_I4();
            LIBS.rotateY(m, angleY);
            LIBS.rotateX(m, angleX);
            LIBS.set_position(m, x, y, z_final);
            return m;
        };

        const headTexture = this.renderer.loadTexture("Resource/prinplup_texture_bare.png");
        const handTexture = this.renderer.loadTexture("Resource/prinplup_hand_texture.png");

        const body_profile = [
            [0.0, 1.95, 0], [0.15, 1.95, 0], [0.3, 1.9, 0], [0.4, 1.8, 0], [0.475, 1.7, 0],
            [0.5, 1.6, 0], [0.525, 1.5, 0], [0.53, 1.4, 0], [0.535, 1.3, 0], [0.55, 1.2, 0],
            [0.57, 1.1, 0], [0.59, 1.0, 0], [0.62, 0.9, 0], [0.65, 0.8, 0], [0.70, 0.7, 0],
            [0.75, 0.6, 0], [0.80, 0.5, 0], [0.82, 0.4, 0], [0.84, 0.3, 0], [0.86, 0.2, 0],
            [0.88, 0.1, 0], [0.92, 0, 0], [0.94, -0.1, 0], [0.96, -0.2, 0], [0.96, -0.3, 0],
            [0.96, -0.4, 0], [0.96, -0.5, 0], [0.96, -0.6, 0], [0.96, -0.7, 0], [0.96, -0.8, 0],
            [0.96, -0.9, 0], [0.90, -1.0, 0], [0.84, -1.1, 0], [0.64, -1.2, 0], [0.44, -1.21, 0],
            [0.24, -1.22, 0], [0.04, -1.23, 0], [0.0, -1.24, 0],
        ];

        // --- Build the Hierarchy ---

        // 1. Create the invisible "breathing" node.
        // All parts that move up/down together will be a child of this.
        const breathingNode = new ModelNode(gl);
        this.rootNode.addChild(breathingNode);
        this.animatedNodes.breathingNode = breathingNode;

        // 2. Create and add all parts as children of 'breathingNode'

        // Body
        const bodyNode = new ModelNode(gl, Geometry.generateLathe(body_profile, 30, C.BODY), headTexture);
        const bodyTrans = (() => {
            const m = createTransform(0, 0.2, 0);
            LIBS.rotateY(m, Math.PI / 2);
            LIBS.scale(m, 1.2);
            return m
        })();
        bodyNode.setBaseTransform(bodyTrans);
        breathingNode.addChild(bodyNode);

        // Body Decorations
        const deco1 = new ModelNode(gl, Geometry.generateCircle(0.17, 20, C.WHITE));
        deco1.setBaseTransform(createOrientedTransform(0.6, 1.1, 1.15, -0.4, 0.25));
        breathingNode.addChild(deco1);

        const deco2 = new ModelNode(gl, Geometry.generateCircle(0.17, 20, C.WHITE));
        deco2.setBaseTransform(createOrientedTransform(0.6, 1.1, 1.15, 0.4, 0.25));
        breathingNode.addChild(deco2);

        const deco3 = new ModelNode(gl, Geometry.generateCircle(0.17, 20, C.WHITE));
        deco3.setBaseTransform(createTransform(-0.4, -0.4, 0.94));
        breathingNode.addChild(deco3);

        const deco4 = new ModelNode(gl, Geometry.generateCircle(0.17, 20, C.WHITE));
        deco4.setBaseTransform(createTransform(0.4, -0.4, 0.94));
        breathingNode.addChild(deco4);

        // Head Disks
        const disk1 = new ModelNode(gl, Geometry.generateSphere(0.7, 0.12, 0.55, 20, 20, C.BEAK));
        disk1.setBaseTransform((() => {
            let m = createTransform(-0.2, 2.1, 0);
            LIBS.rotateZ(m, LIBS.degToRad(120));
            return m;
        })());
        breathingNode.addChild(disk1);

        const disk2 = new ModelNode(gl, Geometry.generateSphere(0.7, 0.12, 0.55, 20, 20, C.BEAK));
        disk2.setBaseTransform((() => {
            let m = createTransform(0.2, 2.1, 0);
            LIBS.rotateZ(m, LIBS.degToRad(60));
            return m;
        })());
        breathingNode.addChild(disk2);

        // Eyes (These are animated, so we save them)
        const leftEye1 = new ModelNode(gl, Geometry.generateAngryEye(0.15, 0.2, 0.1, 10, 10, 130, C.WHITE));
        leftEye1.setBaseTransform((() => {
            let m = createTransform(-0.3, 1.9, 0.4);
            LIBS.rotateZ(m, LIBS.degToRad(330));
            return m;
        }) ());
        breathingNode.addChild(leftEye1);
        this.animatedNodes.eyeNodes.push(leftEye1);

        const leftEye2 = new ModelNode(gl, Geometry.generateAngryEye(0.1, 0.15, 0.12, 10, 10, 130, C.HEAD));
        leftEye2.setBaseTransform((() => {
            let m = createTransform(-0.33, 1.9, 0.42);
            LIBS.rotateZ(m, LIBS.degToRad(330));
            return m;
        }) ());
        breathingNode.addChild(leftEye2);
        this.animatedNodes.eyeNodes.push(leftEye2);

        const leftEye3 = new ModelNode(gl, Geometry.generateAngryEye(0.07, 0.12, 0.12, 10, 10, 130, C.BLACK));
        leftEye3.setBaseTransform((() => {
            let m = createTransform(-0.3, 1.9, 0.45);
            LIBS.rotateZ(m, LIBS.degToRad(330));
            return m;
        }) ());
        breathingNode.addChild(leftEye3);
        this.animatedNodes.eyeNodes.push(leftEye3);

        const leftEye4 = new ModelNode(gl, Geometry.generateSphere(0.03, 0.03, 0.07, 10, 10, C.WHITE));
        leftEye4.setBaseTransform(createTransform(-0.3, 1.93, 0.52));
        breathingNode.addChild(leftEye4);
        this.animatedNodes.eyeNodes.push(leftEye4);

        // Right Eye
        const rightEye1 = new ModelNode(gl, Geometry.generateAngryEye(0.15, 0.2, 0.1, 10, 10, 130, C.WHITE));
        rightEye1.setBaseTransform((() => {
            let m = createTransform(0.3, 1.9, 0.4);
            LIBS.rotateZ(m, LIBS.degToRad(30));
            return m;
        }) ());
        breathingNode.addChild(rightEye1);
        this.animatedNodes.eyeNodes.push(rightEye1);

        const rightEye2 = new ModelNode(gl, Geometry.generateAngryEye(0.1, 0.15, 0.12, 10, 10, 130, C.HEAD));
        rightEye2.setBaseTransform((() => {
            let m = createTransform(0.33, 1.9, 0.42);
            LIBS.rotateZ(m, LIBS.degToRad(30));
            return m;
        }) ());
        breathingNode.addChild(rightEye2);
        this.animatedNodes.eyeNodes.push(rightEye2);

        const rightEye3 = new ModelNode(gl, Geometry.generateAngryEye(0.07, 0.12, 0.12, 10, 10, 130, C.BLACK));
        rightEye3.setBaseTransform((() => {
            let m = createTransform(0.3, 1.9, 0.45);
            LIBS.rotateZ(m, LIBS.degToRad(30));
            return m;
        }) ());
        breathingNode.addChild(rightEye3);
        this.animatedNodes.eyeNodes.push(rightEye3);

        const rightEye4 = new ModelNode(gl, Geometry.generateSphere(0.03, 0.03, 0.07, 10, 10, C.WHITE));
        rightEye4.setBaseTransform(createTransform(0.3, 1.93, 0.52));
        breathingNode.addChild(rightEye4);
        this.animatedNodes.eyeNodes.push(rightEye4);

        // Beak
        const beak = new ModelNode(gl, Geometry.generateBeak(0.22, 0.35, 0.55, 20, C.BEAK));
        beak.setBaseTransform((() => {
            let m = createTransform(0, 1.8, 0.9);
            LIBS.rotateX(m, LIBS.degToRad(5));
            return m;
        })());
        breathingNode.addChild(beak);

        const beakCone = new ModelNode(gl, Geometry.generateCone(0.15, 0.3, 15, C.BEAK));
        beakCone.setBaseTransform((() => {
            let m = createTransform(0, 1.95, 0.6);
            LIBS.rotateX(m, LIBS.degToRad(-15));
            return m;
        })());
        breathingNode.addChild(beakCone);

        const breathEffectNode = new ModelNode(gl, Geometry.generateSphere(0.02, 0.02, 0.03, 10, 10, C.WHITE));
        // This transform is RELATIVE TO THE BEAK NODE
        breathEffectNode.setBaseTransform(createTransform(0, 0.04, 0.06));
        breathEffectNode.alpha = 0.0; // Start invisible
        beak.addChild(breathEffectNode); // Attach to the main beak part
        this.animatedNodes.breathEffects.push(breathEffectNode);
        const breathEffectNode2 = new ModelNode(gl, Geometry.generateSphere(0.03, 0.03, 0.04, 10, 10, C.WHITE));
        // This transform is RELATIVE TO THE BEAK NODE
        breathEffectNode2.setBaseTransform(createTransform(0, 0.08, 0.1));
        breathEffectNode2.alpha = 0.0; // Start invisible
        beak.addChild(breathEffectNode2); // Attach to the main beak part
        this.animatedNodes.breathEffects.push(breathEffectNode2);

        // Hands
        const leftHand = new ModelNode(gl, Geometry.generateSphere(0.2, 1.5, 0.5, 15, 15, C.TAIL), handTexture);
        // NEW: Store the matrix and the node
        const leftHandMatrix = (() => {
            let m = createTransform(-1.2, 0.42, 0.1);
            LIBS.rotateZ(m, LIBS.degToRad(-40));
            LIBS.rotateX(m, LIBS.degToRad(-10));
            return m;
        })();
        leftHand.setBaseTransform(leftHandMatrix);
        this.baseTransforms.leftHand = leftHandMatrix;
        this.animatedNodes.leftHand = leftHand;
        breathingNode.addChild(leftHand);

        const rightHand = new ModelNode(gl, Geometry.generateSphere(0.2, 1.5, 0.5, 15, 15, C.TAIL), handTexture);
        // NEW: Store the matrix and the node
        const rightHandMatrix = (() => {
            let m = createTransform(1.2, 0.42, 0.1);
            LIBS.rotateZ(m, LIBS.degToRad(40));
            LIBS.rotateX(m, LIBS.degToRad(-10));
            return m;
        })();
        rightHand.setBaseTransform(rightHandMatrix);
        this.baseTransforms.rightHand = rightHandMatrix;
        this.animatedNodes.rightHand = rightHand;
        breathingNode.addChild(rightHand);

        // === GRUP KAKI KIRI ===
        const leftLegGroup = new ModelNode(gl); 
        const leftLegGroupBaseMatrix = createTransform(-0.65, -1.2, 0.1); // Posisi global grup kaki
        leftLegGroup.setBaseTransform(leftLegGroupBaseMatrix);
        this.rootNode.addChild(leftLegGroup); // <<< TAMBAHKAN KE ROOTNODE (bukan breathingNode)
        this.animatedNodes.leftLegGroup = leftLegGroup; // Simpan referensi
        this.baseTransforms.leftLegGroup = leftLegGroupBaseMatrix; // Simpan base matrix

        // Tungkai Kiri (Anak dari leftLegGroup)
        const leftLeg = new ModelNode(gl, Geometry.generateSphere(0.3, 0.8, 0.4, 10, 10, C.BODY));
        leftLeg.setBaseTransform(createTransform(0, 2.6, 0)); // Posisi 0,0,0 relatif ke grup
        leftLegGroup.addChild(leftLeg); 

        // Telapak Kaki Kiri 1 (Anak dari leftLegGroup)
        // Posisi LOKAL relatif ke grup kaki (-0.65, -1.2, 0.1)
        // Global: (-0.65, -1.9, 0.24) -> Lokal: (0, -0.7, 0.14)
        const leftFoot1 = new ModelNode(gl, Geometry.generateSphere(0.3, 0.12, 0.35, 10, 10, C.FEET));
        leftFoot1.setBaseTransform(createTransform(0, 1.9, 0.1)); 
        leftLegGroup.addChild(leftFoot1);

        // Telapak Kaki Kiri 2 (Anak dari leftLegGroup)
        // Global: (-0.65, -2.0, 0.4) -> Lokal: (0, -0.8, 0.3)
        const leftFoot2 = new ModelNode(gl, Geometry.generateSphere(0.3, 0.1, 0.5, 10, 10, C.FEET));
        leftFoot2.setBaseTransform(createTransform(0, 1.8, 0.25));
        leftLegGroup.addChild(leftFoot2);

        // === GRUP KAKI KANAN ===
        const rightLegGroup = new ModelNode(gl);
        const rightLegGroupBaseMatrix = createTransform(0.65, -1.2, 0.1); // Posisi global grup kaki
        rightLegGroup.setBaseTransform(rightLegGroupBaseMatrix);
        this.rootNode.addChild(rightLegGroup); // <<< TAMBAHKAN KE ROOTNODE
        this.animatedNodes.rightLegGroup = rightLegGroup;
        this.baseTransforms.rightLegGroup = rightLegGroupBaseMatrix;

        // Tungkai Kanan (Anak dari rightLegGroup)
        const rightLeg = new ModelNode(gl, Geometry.generateSphere(0.3, 0.8, 0.4, 10, 10, C.BODY));
        rightLeg.setBaseTransform(createTransform(0, 2.6, 0)); 
        rightLegGroup.addChild(rightLeg);

        // Telapak Kaki Kanan 1 (Anak dari rightLegGroup)
        // Global: (0.65, -1.9, 0.24) -> Lokal: (0, -0.7, 0.14)
        const rightFoot1 = new ModelNode(gl, Geometry.generateSphere(0.3, 0.12, 0.35, 10, 10, C.FEET));
        rightFoot1.setBaseTransform(createTransform(0, 1.9, 0.1));
        rightLegGroup.addChild(rightFoot1);

        // Telapak Kaki Kanan 2 (Anak dari rightLegGroup)
        // Global: (0.65, -2.0, 0.4) -> Lokal: (0, -0.8, 0.3)
        const rightFoot2 = new ModelNode(gl, Geometry.generateSphere(0.3, 0.1, 0.6, 10, 10, C.FEET));
        rightFoot2.setBaseTransform(createTransform(0, 1.8, 0.25));
        rightLegGroup.addChild(rightFoot2);

        // Tail
        const tail = new ModelNode(gl, Geometry.generateTaperedShapeFromSpline(
            [[0.0, -0, -0.2], [0.0, -0.9, -1.2], [0.0, -1.2, -2.0]],
            50, [0.9, 0.6], [0.01, 0.01], 20, C.TAIL
        ));
        tail.setBaseTransform(createTransform(0, -0.2, 0.1));
        breathingNode.addChild(tail);
    }

    // NEW: Central animation update function
    updateAnimation(animValues) {
        // 1. Apply 'bodyBreathe' Y-translation to the main breathing node
        const T_breath = LIBS.get_I4();
        LIBS.translateY(T_breath, animValues.bodyTranslate + 2.5);
        LIBS.scaleY(T_breath, animValues.bodyScale);
        LIBS.scale(T_breath, animValues.bodyScale - 0.01);
        this.animatedNodes.breathingNode.setLocalTransform(T_breath);

        // 2. Apply ONLY 'bodyTranslate' to the leg groups
        const T_legBreath = LIBS.get_I4();
        LIBS.translateY(T_legBreath, animValues.bodyTranslate || 0.0); // Hanya translasi

        // Terapkan ke Kaki Kiri
        const leftLegNode = this.animatedNodes.leftLegGroup;
        if (leftLegNode) {
            // M_final = M_base * M_animasi_translasi
            const finalLeftLegMatrix = LIBS.multiply(this.baseTransforms.leftLegGroup, T_legBreath);
            leftLegNode.setLocalTransform(finalLeftLegMatrix);
        }

        // Terapkan ke Kaki Kanan
        const rightLegNode = this.animatedNodes.rightLegGroup;
        if (rightLegNode) {
            // M_final = M_base * M_animasi_translasi
            const finalRightLegMatrix = LIBS.multiply(this.baseTransforms.rightLegGroup, T_legBreath);
            rightLegNode.setLocalTransform(finalRightLegMatrix);
        }

        // (Future animations like 'diskBreathe' would be added here)

        // 3. NEW: Apply 'flapAngle' Z-rotation to hands
        const flapAngle = animValues.flapAngle || 0.0;
        const pivotY = 1.5; // From the hand's geometry radius
        const pivotX = 0.7;

        // Create pivot matrices
        let T_up = LIBS.get_I4();
        LIBS.translateY(T_up, -pivotY);
        LIBS.translateX(T_up, pivotX);
        let T_down = LIBS.get_I4();
        LIBS.translateY(T_down, pivotY);
        LIBS.translateX(T_down, -pivotX);

        // --- Left Hand ---
        const R_left = LIBS.get_I4();
        LIBS.rotateZ(R_left, flapAngle);

        // Combine for animation matrix: M_anim = T_up * R_z * T_down
        let leftAnim = LIBS.multiply(R_left, T_down);
        leftAnim = LIBS.multiply(T_up, leftAnim);

        // Combine with base: M_final = M_base * M_anim
        const leftFinal = LIBS.multiply(this.baseTransforms.leftHand, leftAnim);
        this.animatedNodes.leftHand.setLocalTransform(leftFinal);

        // --- Right Hand ---
        // Create pivot matrices
        T_up = LIBS.get_I4();
        LIBS.translateY(T_up, -pivotY);
        LIBS.translateX(T_up, -pivotX);
        T_down = LIBS.get_I4();
        LIBS.translateY(T_down, pivotY);
        LIBS.translateX(T_down, pivotX);

        const R_right = LIBS.get_I4();
        LIBS.rotateZ(R_right, -flapAngle); // Opposite direction

        // Combine for animation matrix: M_anim = T_up * R_z * T_down
        let rightAnim = LIBS.multiply(R_right, T_down);
        rightAnim = LIBS.multiply(T_up, rightAnim);

        // Combine with base: M_final = M_base * M_anim
        const rightFinal = LIBS.multiply(this.baseTransforms.rightHand, rightAnim);
        this.animatedNodes.rightHand.setLocalTransform(rightFinal);

        const breathAlpha = animValues.breathAlpha || 0.0;
        const breathScale = animValues.breathScale || 1.0;

        // Buat matriks skala sekali saja
        const S_breath = LIBS.get_I4();
        LIBS.scale(S_breath, breathScale);

        // Loop melalui setiap node efek napas di array
        this.animatedNodes.breathEffects.forEach(breathNode => {
            if (breathNode) {
                // Update Alpha
                breathNode.alpha = breathAlpha;

                // Gabungkan base transform node DENGAN matriks skala
                // M_local = M_base * M_animasi_skala
                const finalBreathMatrix = LIBS.multiply(breathNode.baseMatrix, S_breath);
                breathNode.setLocalTransform(finalBreathMatrix);
            }
        });
    }

    draw(shader, parentMatrix) {
        // 'parentMatrix' is the animated matrix from the ice island
        // 'this.modelMatrix' is the mouse-drag rotation
        const finalParentMatrix = LIBS.multiply(parentMatrix, this.modelMatrix);

        // Update all world matrices starting from the root
        this.rootNode.updateWorldMatrix(finalParentMatrix);

        // Start the recursive draw
        this.rootNode.draw(shader);
    }
}