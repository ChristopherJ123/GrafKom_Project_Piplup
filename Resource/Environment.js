import { Geometry } from "./Geometries.js";
import { LIBS } from "./Libs.js";
import { ModelNode } from "./ModelNode.js";

export class Environment {
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

        const waterTexture = this.renderer.loadTexture('Resource/water-texture.png');

        // 1. Ice Island
        const iceIslandNode = new ModelNode(gl, Geometry.generateIrregularExtrudedPolygon(8, 15, 1, C.SNOW_WHITE, 0.5));
        iceIslandNode.setBaseTransform(LIBS.get_I4());
        this.rootNode.addChild(iceIslandNode);
        this.animatedNodes.iceIsland = iceIslandNode; // Save for animation

        // 2. Water
        const waterNode = new ModelNode(gl, Geometry.generateWaterPlane(200, 200, 1, 1, C.WATER_BLUE), waterTexture);
        waterNode.setBaseTransform(LIBS.get_I4());
        this.rootNode.addChild(waterNode);
    }

    // NEW: Central animation update
    updateAnimation() {
        this.animationTime += 0.02;

        // 1. Calculate float animation
        const amplitude = 0.1;
        const floatY = Math.sin(this.animationTime) * amplitude;

        // 2. Apply animation to the ice island
        // Note: We combine with the base matrix
        const T_float = LIBS.get_I4();
        LIBS.translateY(T_float, floatY);

        const iceNode = this.animatedNodes.iceIsland;
        const finalIceMatrix = LIBS.multiply(iceNode.baseMatrix, T_float);
        iceNode.setLocalTransform(finalIceMatrix);
    }

    // MODIFIED: Return the final computed world matrix of the island
    getIceIslandWorldMatrix() {
        // Ensure the matrix is up-to-date before returning it
        this.rootNode.updateWorldMatrix(this.modelMatrix);
        return this.animatedNodes.iceIsland.worldMatrix;
    }

    draw(shader) {
        // Update all world matrices
        this.rootNode.updateWorldMatrix(this.modelMatrix);
        // Start recursive draw
        this.rootNode.draw(shader);
    }
}