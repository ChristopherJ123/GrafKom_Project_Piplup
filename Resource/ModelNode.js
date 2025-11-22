import { LIBS } from "./Libs.js";

export class ModelNode {
    constructor(gl, geometry = null, texture = null) {
        this.gl = gl;
        this.geometry = geometry; // The drawable geometry
        this.texture = texture;
        this.buffers = null;
        this.alpha = 1.0;

        this.baseMatrix = LIBS.get_I4();  // The "default" pose (set once)
        this.localMatrix = LIBS.get_I4(); // The final transform (base * animation)
        this.worldMatrix = LIBS.get_I4(); // Final transformation in world space
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

    createBuffers() {
        const vertexBuffer = this.gl.createBuffer();
        this.gl.bindBuffer(this.gl.ARRAY_BUFFER, vertexBuffer);
        this.gl.bufferData(this.gl.ARRAY_BUFFER, new Float32Array(this.geometry.vertices), this.gl.STATIC_DRAW);

        const facesBuffer = this.gl.createBuffer();
        this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, facesBuffer);
        this.gl.bufferData(this.gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(this.geometry.faces), this.gl.STATIC_DRAW);

        return { vertex: vertexBuffer, faces: facesBuffer, faces_length: this.geometry.faces.length };
    }

    // Set this node's *base* transform (e.g., the flipper's position relative to the body)
    setBaseTransform(matrix) {
        this.baseMatrix = matrix;
        this.localMatrix = matrix; // By default, local = base
    }

    // Set this node's *animated* local transform
    setLocalTransform(matrix) {
        this.localMatrix = matrix;
    }

    // Recursive function to update all matrices in the tree
    updateWorldMatrix(parentWorldMatrix) {
        if (parentWorldMatrix) {
            this.worldMatrix = LIBS.multiply(this.localMatrix, parentWorldMatrix);
        } else {
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
            gl.uniformMatrix4fv(shader.locations.Mmatrix, false, this.worldMatrix);
            gl.uniform1f(shader.locations.u_alpha, this.alpha);

            gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.vertex);

            // Stride is now 11 floats (3 pos, 3 color, 2 uv, 3 normal)
            const stride = 4 * (3 + 3 + 2 + 3);
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

        // Now, recursively draw all children
        for (const child of this.children) {
            child.draw(shader);
        }
    }
}