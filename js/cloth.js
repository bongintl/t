function Cloth(src, x, y, z, w, h){
	
	var cloth = this,
		restDistance = 1,
		damping = 0.5,
		restDistanceSq = restDistance * restDistance,
		windMultiplier = .1
		xSegs = Math.floor(w / restDistance),
		ySegs = Math.floor(h / restDistance);
	
	cloth.src = src;
	cloth.w = w - (w % restDistance);
	cloth.h = h - (h % restDistance);
	cloth.x = x;
	cloth.y = y - cloth.h;
	cloth.z = z;
	
	function geometryFunction(u, v){
		return new THREE.Vector3( cloth.x + u * cloth.w, cloth.y + v * cloth.h, cloth.z );
	}
	
	cloth.geometry = new THREE.ParametricGeometry(geometryFunction, xSegs, ySegs);
	cloth.geometry.dynamic = true;
	
	var texture = THREE.ImageUtils.loadTexture( src );
	texture.minFilter = THREE.LinearFilter;
	
	cloth.material = new THREE.MeshLambertMaterial({
		//wireframe: true,
        map: texture,
        side: THREE.DoubleSide
    });
    
    cloth.mesh = new THREE.Mesh(this.geometry, this.material);
    cloth.vertices = cloth.geometry.vertices;
	cloth.vertices2 = [];
    cloth.forces = [];
    cloth.vertexFaces = [];
	cloth.pinnedVertices = [];
    
    cloth.arrayXY = function (array, x, y){
	    return array[(xSegs+1) * y + x];
    }
        
    cloth.getVertex = function(x, y){
	    return cloth.arrayXY(cloth.geometry.vertices, x, y);
    }
    
    for(var i = 0, ii = cloth.vertices.length; i < ii; ++i){
	    cloth.forces[i] = new THREE.Vector3();
	    cloth.vertexFaces[i] = 0;
		cloth.vertices2.push(cloth.vertices[i].clone());
		cloth.pinnedVertices[i] = i >= (cloth.vertices.length - xSegs - 1);
    }
    	
	var face;
	for(i = 0, ii = cloth.geometry.faces.length; i < ii; ++i){
		
		face = cloth.geometry.faces[i];
		
		cloth.vertexFaces[face.a]++;
		cloth.vertexFaces[face.b]++;
		cloth.vertexFaces[face.c]++;
		
	}
		
	cloth.move = function(v){
				
		var pinnedVertex;
		
		cloth.x += v.x;
		cloth.y += v.y;
		cloth.z += v.z;
		
		cloth.mesh.position.add(v);
		
		for(var i = 0, ii = xSegs; i <= ii; ++i){
			pinnedVertex = cloth.getVertex(i, ySegs);
			pinnedVertex.add(v);
		}
		
	}
	
	cloth.moveTo = function(v){
				
		var d = Cloth.tmpVec3.copy(v).sub(cloth.mesh.position);
		cloth.move(d);
		
	}
	
	cloth.integrate = function(frame){
				
		var v1, v2, f;
		
		
	}
		
	cloth.tick = function(frustum, gravity, wind){
									
		var	onScreen,
			i, ii,
			j, jj,
			face,
			normal,
			vertex,
			v1, v2, v3,
			appliedWindForce,
			delta,
			dist,
			distSq,
			correction,
			f;
			
		for(i = 0, ii = cloth.vertices.length; i < ii; ++i){
			if(frustum.containsPoint(cloth.vertices[i])) {
				onScreen = true;
				break;
			}
		}
		
		if(!onScreen) return;
				
		// Add gravity to unpinned vertices	
		for(i = 0, ii = cloth.vertices.length; i < ii; ++i){
			
			if(cloth.pinnedVertices[i]) continue;
			cloth.forces[i].add(gravity);

		}
		
		if(!onScreen) return;
		
		// Find wind force on faces and add to unpinned vertices
		for(i = 0, ii = cloth.geometry.faces.length; i < ii; ++i){
						
			face = cloth.geometry.faces[i];
			normal = face.normal;
			appliedWindForce = Cloth.tmpVec3
				.copy(normal)
				.multiplyScalar(normal.dot(wind));
				
			v1 = cloth.forces[face.a];
			v2 = cloth.forces[face.b];
			v3 = cloth.forces[face.c];
			
			if(!cloth.pinnedVertices[face.a]) v1.add(appliedWindForce.clone().multiplyScalar(windMultiplier/cloth.vertexFaces[face.a]));
			if(!cloth.pinnedVertices[face.b]) v2.add(appliedWindForce.clone().multiplyScalar(windMultiplier/cloth.vertexFaces[face.b]));
			if(!cloth.pinnedVertices[face.c]) v3.add(appliedWindForce.clone().multiplyScalar(windMultiplier/cloth.vertexFaces[face.c]));
			
		}
		
			
		// Satisfy vertical constraints
		for(i = 0, ii = xSegs; i <= ii; ++i){
			
			for(j = ySegs, jj = 0; j > jj; --j){
								
				v1 = cloth.getVertex(i, j);
				v2 = cloth.getVertex(i, j-1);
				
				delta = Cloth.tmpVec3.subVectors( v2, v1 );
				distSq = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
												
				correction = Cloth.tmpVec3.multiplyScalar( restDistanceSq / ( distSq + restDistanceSq ) -.5);
				
				v2.add(correction);
				
			}
			
		}

		// Satisfy horizontal constraints
		for(i = 1, ii = ySegs; i <= ii; ++i){
			
			for(j = 0, jj = xSegs; j < jj; ++j){
								
				v1 = cloth.getVertex(j, i);
				v2 = cloth.getVertex(j+1, i);
				
				delta = Cloth.tmpVec3.subVectors( v2, v1 );
				distSq = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
								
				correction = Cloth.tmpVec3.multiplyScalar( restDistanceSq / ( distSq + restDistanceSq ) - .5 );
				
				v1.sub(correction);
				v2.add(correction);
				
			}
			
		}
		
		// Verlet integration
		for(i = 0, ii = cloth.vertices.length; i < ii; ++i){
			
			if(cloth.pinnedVertices[i]) continue;
			
			v1 = cloth.vertices[i];
			v2 = cloth.vertices2[i];
			f = cloth.forces[i];
						
			Cloth.tmpVec3.copy(v1).multiplyScalar(2-damping).sub(v2.clone().multiplyScalar(1-damping)).add(f.multiplyScalar(1));
			
			v2.copy(v1);
			v1.copy(Cloth.tmpVec3);
			
			f.set(0,0,0);
			
		}
								
		cloth.geometry.computeFaceNormals();
		cloth.geometry.computeVertexNormals();
		cloth.geometry.computeBoundingBox();
		
		cloth.geometry.normalsNeedUpdate = true;
		cloth.geometry.verticesNeedUpdate = true;
		
	}
	
};

Cloth.tmpVec3 = new THREE.Vector3();