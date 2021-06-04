# Copyright (C) 2010-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import espressomd
import numpy as np
import ipdb

class Mesh:
    '''
    mesh data of elastic objects
    '''
    def __init__(self):
        # input files
        self.nodes_file = ""
        self.triangles_file = ""

        # mesh data
        self.vertices = []
        self.edges = []
        self.faces = []
        self.dihedrals = []

        self.vertices_num = 0
        self.edges_num = 0
        self.faces_num = 0
        self.dihedrals_num = 0

        # look-up dicts
        self.edges_vertices2index = {}
        self.faces_vertices2index = {}
        self.dihedrals_vertices2index = {}

    def load_nodes(self, filename):
        '''
        load cooridnates of each vertex
        '''
        self.vertices_num = 0

        with open(filename) as nodesfile:
            for line in nodesfile:
                data = line.split()

                index = int(data[0])
                x = float(data[1])
                y = float(data[2])
                z = float(data[3])

                self.vertices.append(Vertex(index, float(x), float(y), float(z)))

                self.vertices_num += 1

        print('{} vertices loaded from nodes file: {}'.format(self.vertices_num, filename))

    def load_triangles(self, filename):
        '''
        load triangles of the surface
        '''
        self.faces_num = 0

        with open(filename) as trianglesfile:
            for line in trianglesfile:
                data = line.split()

                index = int(data[0])
                v1 = int(data[1])
                v2 = int(data[2])
                v3 = int(data[3])

                r1 = self.vertices[v1].pos
                r2 = self.vertices[v2].pos
                r3 = self.vertices[v3].pos

                # assume nodes are centered around origin
                ort = int(np.sign( np.dot(r1, np.cross(r2, r3)) ))

                if ort == 1:
                    newFace = Face(index, v1, v2, v3)
                    self.faces.append(newFace)
                    self.faces_vertices2index.update({newFace.vertices:index})
                else:
                    newFace = Face(index, v3, v2, v1)
                    self.faces.append(newFace)
                    self.faces_vertices2index.update({newFace.vertices:index})

                self.faces_num += 1

        print('{} faces loaded from triangles file: {}'.format(self.faces_num, filename))

    def build_edges(self):
        '''
        build edges according to faces
        '''
        edges = set()

        for face in self.faces:
            edge1 = frozenset([face.v1, face.v2])
            edge2 = frozenset([face.v2, face.v3])
            edge3 = frozenset([face.v3, face.v1])

            edges.add(edge1)
            edges.add(edge2)
            edges.add(edge3)

        self.edges_num = 0

        for edge in edges:
            index = self.edges_num
            vertex_pair = [vertex for vertex in edge]
            v1 = vertex_pair[0]
            v2 = vertex_pair[1]

            newEdge = Edge(index, v1, v2)
            self.edges.append(newEdge)
            self.edges_vertices2index.update({newEdge.vertices:index})

            self.edges_num += 1

        print('{} edges added.'.format(self.edges_num))

    def build_dihedrals(self):
        '''
        build dihedrals according to edges
        '''
        self.dihedrals_num = 0

        for edge in self.edges:
            edge_vertices = edge.vertices

            # find two faces shared this edge
            adjacent_faces = []
            for face in self.faces:
                if edge_vertices.issubset(face.vertices):
                    adjacent_faces.append(face)
            assert len(adjacent_faces) == 2

            face1 = adjacent_faces[0]
            face2 = adjacent_faces[1]

            index = self.dihedrals_num
            v1 = next(iter(face1.vertices.difference(edge_vertices)))
            v2 = edge.v1
            v3 = edge.v2
            v4 = next(iter(face2.vertices.difference(edge_vertices)))

            newDihedral = Dihedral(index, v1, v2, v3, v4)
            self.dihedrals.append(newDihedral)
            self.dihedrals_vertices2index.update({newDihedral.vertices:index})

            self.dihedrals_num += 1

        print('{} dihedrals added.'.format(self.dihedrals_num))

class Vertex:
    '''
    represent a vertex data (0D)
    '''
    def __init__(self, index, x, y, z):
        self.index = index
        self.x = x
        self.y = y
        self.z = z
        self.pos = np.array([x,y,z])

        # neighbors info
        self.edges_num = 0
        self.neighbor_vertices = []

class Edge:
    '''
    represent an edge data (1D)
    '''
    def __init__(self, index, v1, v2):
        self.index = index
        self.v1 = v1
        self.v2 = v2
        self.vertices = frozenset([v1, v2])

class Face:
    '''
    represent a face data (2D)
    '''
    def __init__(self, index, v1, v2, v3):
        self.index = index
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.vertices = frozenset([v1, v2, v3])

class Dihedral:
    '''
    represent a dihedral angle
    '''
    def __init__(self, index, v1, v2, v3, v4):
        self.index = index
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.v4 = v4
        self.vertices = frozenset([v1, v2, v3, v4])

if __name__ == "__main__":
    # test purpose
    nodesfile = '1082.6.6.nodes'
    trianglesfile = '1082.6.6.triangles'
    testMesh = Mesh()
    testMesh.load_nodes(nodesfile)
    testMesh.load_triangles(trianglesfile)
    testMesh.build_edges()
    testMesh.build_dihedrals()
    ipdb.set_trace()
