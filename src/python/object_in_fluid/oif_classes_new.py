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
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import ipdb
from scipy.spatial.qhull import Voronoi

class Mesh:
    '''
    mesh data of elastic objects
    '''
    def __init__(self, nodes_file="", triangles_file="", rescale=(1.0, 1.0, 1.0)):
        # input files
        self.nodes_file = nodes_file
        self.triangles_file = triangles_file

        # scaling factors
        self.rescale = rescale

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

        # mesh properties
        self.avg_edge_length = None
        self.total_area = None
        self.total_area_voronoi = None
        self.total_volume = None

    def initialize(self):
        '''
        initialize mesh data
        '''

        # load input data

        if self.nodes_file == "":
            raise Exception("Missing nodes file. Abort.")
        self.load_nodes(self.nodes_file)

        if self.triangles_file == "":
            # construct triangulation according to input nodes
            self.build_faces()
        else:
            # load provided triangulation
            self.load_triangles(self.triangles_file)

        # build edges
        self.build_edges()

        # build dihedrals
        self.build_dihedrals()

        # build vertices neighbor info
        self.build_vertices_neighbor_info()

        # type assignment
        self.assign_vertices_type()
        self.assign_edges_type()
        self.assign_dihedrals_type()

        # additional mesh properties
        self.compute_mesh_properties()

        # output mesh statistics
        self.print_mesh_info()

        #self.plot_points()

    def assign_vertices_type(self):
        '''
        assign A/B components
        '''
        cntA = 0
        cntB = 0

        for vertex in self.vertices:
            if vertex.x > 0.0:
                vertex.set_type(0)
                cntA += 1
            elif vertex.x < 0.0:
                vertex.set_type(1)
                cntB += 1
            else:
                # random for vertex in the middle
                type_rnd = np.random.randint(2)
                vertex.set_type(type_rnd)
                if type_rnd == 0:
                    cntA += 1
                else:
                    cntB += 1

        print('A type vertices number: {}'.format(cntA))
        print('B type vertices number: {}'.format(cntB))

    def assign_edges_type(self):
        '''
        assign edge type according to its vertices
        '''
        cntAA = 0
        cntBB = 0
        cntAB = 0

        for edge in self.edges:
            v1 = edge.v1
            v2 = edge.v2

            v1_type = self.vertices[v1].type
            v2_type = self.vertices[v2].type
            if v1_type == v2_type:
                edge.set_type(v1_type)
                if v1_type == 0:
                    cntAA += 1
                else:
                    cntBB += 1
            else:
                edge.set_type(2)
                cntAB += 1

        print('AA type edges number: {}'.format(cntAA))
        print('BB type edges number: {}'.format(cntBB))
        print('AB type edges number: {}'.format(cntAB))

    def assign_dihedrals_type(self):
        '''
        assign dihedral type according to common edge vertices
        '''
        cntAA = 0
        cntBB = 0
        cntAB = 0

        for dihedral in self.dihedrals:
            v2 = dihedral.v2
            v3 = dihedral.v3

            v2_type = self.vertices[v2].type
            v3_type = self.vertices[v3].type
            if v2_type == v3_type:
                dihedral.set_type(v2_type)
                if v2_type == 0:
                    cntAA += 1
                else:
                    cntBB += 1
            else:
                dihedral.set_type(2)
                cntAB += 1

        print('AA type dihedrals number: {}'.format(cntAA))
        print('BB type dihedrals number: {}'.format(cntBB))
        print('AB type dihedrals number: {}'.format(cntAB))

    def plot_points(self):
        fig = plt.figure()
        ax = Axes3D(fig)
        #vertices = np.zeros( (self.vertices_num, 3) )
        #for i in range(self.vertices_num):
        #    vertices[i][0] = self.vertices[i].x
        #    vertices[i][1] = self.vertices[i].y
        #    vertices[i][2] = self.vertices[i].z
        #ax.scatter(vertices[:,0], vertices[:,1], vertices[:,2])
        dihedral = self.dihedrals[0]
        print(dihedral.index)
        r1 = dihedral.r1
        r2 = dihedral.r2
        r3 = dihedral.r3
        r4 = dihedral.r4

        ax.scatter(dihedral.r1[0], dihedral.r1[1], dihedral.r1[2])
        ax.scatter(dihedral.r2[0], dihedral.r2[1], dihedral.r2[2])
        ax.scatter(dihedral.r3[0], dihedral.r3[1], dihedral.r3[2])
        ax.scatter(dihedral.r4[0], dihedral.r4[1], dihedral.r4[2])

        ax.text(dihedral.r1[0], dihedral.r1[1], dihedral.r1[2], '1: {}'.format(dihedral.v1) )
        ax.text(dihedral.r2[0], dihedral.r2[1], dihedral.r2[2], '2: {}'.format(dihedral.v2) )
        ax.text(dihedral.r3[0], dihedral.r3[1], dihedral.r3[2], '3: {}'.format(dihedral.v3) )
        ax.text(dihedral.r4[0], dihedral.r4[1], dihedral.r4[2], '4: {}'.format(dihedral.v4) )

        ax.plot([r1[0], r2[0]],[r1[1], r2[1]],[r1[2], r2[2]])
        ax.plot([r2[0], r3[0]],[r2[1], r3[1]],[r2[2], r3[2]], ':')
        ax.plot([r3[0], r1[0]],[r3[1], r1[1]],[r3[2], r1[2]])
        ax.plot([r3[0], r4[0]],[r3[1], r4[1]],[r3[2], r4[2]], '--')
        plt.show()

    def compute_mesh_properties(self):
        '''
        compute additional properties based on underlying mesh data
        '''
        self.compute_face_properties()

        self.compute_vertex_properties()

        self.compute_average_edge_length()

        self.compute_total_area()

        self.compute_total_volume()

    def print_mesh_info(self):
        '''
        output a summary of mesh properties
        '''
        print("###########################")
        print("#      Mesh Info          #")
        print("###########################")

        print("vertices number = {} edges number = {} faces number = {}".format(self.vertices_num, self.edges_num, self.faces_num))
        print("total area = {} total volume = {}".format(self.total_area, self.total_volume))
        print("total area (voronoi) = {}".format(self.total_area_voronoi))
        print("average edge length = {}".format(self.avg_edge_length))
        gb = 0.0
        for vertex in self.vertices:
            gb += vertex.area * vertex.gaussian_curvature
        print("Gauss-Bonnet = {} deviation = {}".format(gb, np.abs(4*np.pi-gb)))
    
    def compute_vertex_properties(self):
        '''
        compute vertex properties
        '''
        # normals
        self.compute_vertex_normals()
        # curvatures
        self.compute_vertex_curvatures()

    def compute_vertex_normals(self):
        '''
        area-weighted normal vector
        '''
        for vertex in self.vertices:
            vertex_normal = np.zeros(3)
            normalization = 0.0
            for neighbor_face_index in vertex.neighbor_faces:
                face_area = self.faces[neighbor_face_index].area
                face_normal = self.faces[neighbor_face_index].normal
                vertex_normal += face_normal * face_area
                normalization += face_area
            vertex_normal /= normalization
            vertex.normal = vertex_normal / np.linalg.norm(vertex_normal)

    def compute_vertex_curvatures(self):
        '''
        compute vertex curvatures
        ref: Gompper and Kroll 1996
        ref: Discrete Differential Geometry Operators for Triangulated 2-Manifolds 2003
        '''
        # loop for all vertices
        for vertex in self.vertices:
            voronoi_area = 0.0
            mean_curvature_vec = np.zeros(3)
            internal_angles = 0.0
            # loop for all neighbor vertices
            for neighbor_vertex_index in vertex.neighbor_vertices:
                neighbor_vertex = self.vertices[neighbor_vertex_index]
                edge_vertices = frozenset([vertex.index, neighbor_vertex_index])
                edge_index = self.edges_vertices2index[edge_vertices]
                edge = self.edges[edge_index]

                ri = vertex.r
                rj = neighbor_vertex.r
                rij = rj - ri
                rij_norm = np.linalg.norm(rij)

                opposite_vertices_index = [v for v in edge.opposite_vertices]
                assert len(opposite_vertices_index) == 2
                v1 = opposite_vertices_index[0]
                v2 = opposite_vertices_index[1]
                r1 = self.vertices[v1].r
                r2 = self.vertices[v2].r

                # mean-curvature
                r1i = ri - r1
                r1i_norm = np.linalg.norm(r1i)
                r1j = rj - r1
                r1j_norm = np.linalg.norm(r1j)
                costheta1 = np.dot(r1i/r1i_norm, r1j/r1j_norm)
                cottheta1 = costheta1 / np.sqrt(1-costheta1*costheta1)

                r2i = ri - r2
                r2i_norm = np.linalg.norm(r2i)
                r2j = rj - r2
                r2j_norm = np.linalg.norm(r2j)
                costheta2 = np.dot(r2i/r2i_norm, r2j/r2j_norm)
                cottheta2 = costheta2 / np.sqrt(1-costheta2*costheta2)

                voronoi_area += (cottheta1 + cottheta2) * rij_norm**2 / 8.0
                mean_curvature_vec += (cottheta1 + cottheta2) * rij / 2.0

                # gaussian-curvature
                ri1 = -r1i
                ri1_norm = r1i_norm
                ri2 = -r2i
                ri2_norm = r2i_norm

                theta1 = np.arccos( np.dot(rij/rij_norm, ri1/ri1_norm) )
                if theta1 > np.pi / 2.0:
                    raise Exception("Obtuse angle exists. Correction not implemented yet.")
                internal_angles += theta1

                theta2 = np.arccos( np.dot(rij/rij_norm, ri2/ri2_norm) )
                if theta2 > np.pi / 2.0:
                    raise Exception("Obtuse angle exists. Correction not implemented yet.")
                internal_angles += theta2
            
            # voroni cell area
            vertex.area = voronoi_area

            # mean curvature
            # sign convention:
            # convex region > 0 concave region < 0
            mean_curvature_vec /= voronoi_area
            mean_curvature_norm = np.linalg.norm(mean_curvature_vec)
            ort = np.sign( np.dot(mean_curvature_vec/mean_curvature_norm, vertex.normal) )
            vertex.mean_curvature = -ort * mean_curvature_norm / 2.0

            # gaussian curvature
            internal_angles /= 2.0 # internal angles double counted
            gaussian_curvature = (2.0 * np.pi - internal_angles) / voronoi_area
            vertex.gaussian_curvature = gaussian_curvature

            # principal curvatures
            delta = 0.0
            mean_curvature_sqr = vertex.mean_curvature**2
            if mean_curvature_sqr > gaussian_curvature:
                delta = mean_curvature_sqr - gaussian_curvature
            else:
                # unphysical case due to numerical errors
                delta = 0.0
                #print("WARNING: square of mean curvature is smaller than gaussian curvature.")
            
            vertex.curvature_1 = vertex.mean_curvature + np.sqrt(delta)
            vertex.curvature_2 = vertex.mean_curvature - np.sqrt(delta)

    def compute_face_properties(self):
        '''
        compute face properties
        '''
        for face in self.faces:
            v1 = face.v1
            v2 = face.v2
            v3 = face.v3

            r1 = self.vertices[v1].r
            r2 = self.vertices[v2].r
            r3 = self.vertices[v3].r

            r12 = r2 - r1
            r13 = r3 - r1

            n = np.cross(r12, r13)
            n_norm = np.linalg.norm(n)

            face.area = n_norm / 2.0
            face.volume = np.dot(r1, np.cross(r2, r3)) / 6.0
            face.normal = n / n_norm

    def compute_average_edge_length(self):
        '''
        average edge length
        '''
        self.avg_edge_length = 0.0
        for edge in self.edges:
            self.avg_edge_length += edge.get_length()
        self.avg_edge_length /= self.edges_num
    
    def compute_total_area(self):
        '''
        total area of mesh
        '''
        self.total_area = 0.0
        for face in self.faces:
            self.total_area += face.area

        self.total_area_voronoi = 0.0
        for vertex in self.vertices:
            self.total_area_voronoi += vertex.area
    
    def compute_total_volume(self):
        '''
        total volume of mesh
        '''
        self.total_volume = 0.0
        for face in self.faces:
            self.total_volume += face.volume

    def load_nodes(self, filename):
        '''
        load cooridnates of each vertex
        '''
        self.vertices_num = 0

        with open(filename) as nodesfile:
            for line in nodesfile:
                data = line.split()

                index = int(data[0])
                x = self.rescale[0] * float(data[1])
                y = self.rescale[1] * float(data[2])
                z = self.rescale[2] * float(data[3])

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

                r1 = self.vertices[v1].r
                r2 = self.vertices[v2].r
                r3 = self.vertices[v3].r

                # assume nodes are centered around origin
                ort = int(np.sign( np.dot(r1, np.cross(r2, r3)) ))

                if ort == 1:
                    newFace = Face(index, v1, v2, v3)
                    newFace.set_vertices_pos(r1, r2, r3)
                    self.faces.append(newFace)
                    self.faces_vertices2index.update({newFace.vertices:index})
                else:
                    newFace = Face(index, v3, v2, v1)
                    newFace.set_vertices_pos(r3, r2, r1)
                    self.faces.append(newFace)
                    self.faces_vertices2index.update({newFace.vertices:index})

                self.faces_num += 1

        print('{} faces loaded from triangles file: {}'.format(self.faces_num, filename))
    
    def build_vertices_neighbor_info(self):
        '''
        get neighbor vertices, edges and faces of each vertex
        '''
        for edge in self.edges:
            index = edge.index
            v1 = edge.v1
            v2 = edge.v2

            # vertex neighbor vertices info
            self.vertices[v1].neighbor_vertices.add(v2)
            self.vertices[v2].neighbor_vertices.add(v1)

            # vertex neighbor edges info
            self.vertices[v1].neighbor_edges.add(index)
            self.vertices[v2].neighbor_edges.add(index)

        for face in self.faces:
            index = face.index
            v1 = face.v1
            v2 = face.v2
            v3 = face.v3

            # vertex neighbor faces info
            self.vertices[v1].neighbor_faces.add(index)
            self.vertices[v2].neighbor_faces.add(index)
            self.vertices[v3].neighbor_faces.add(index)

    def build_edges(self):
        '''
        build edges according to faces
        '''
        edges_set = set()

        for face in self.faces:
            edge1 = frozenset([face.v1, face.v2])
            edge2 = frozenset([face.v2, face.v3])
            edge3 = frozenset([face.v3, face.v1])

            edges_set.add(edge1)
            edges_set.add(edge2)
            edges_set.add(edge3)

        self.edges_num = 0

        for edge in edges_set:
            index = self.edges_num
            vertex_pair = [vertex for vertex in edge]
            v1 = vertex_pair[0]
            v2 = vertex_pair[1]

            r1 = self.vertices[v1].r
            r2 = self.vertices[v2].r

            newEdge = Edge(index, v1, v2)
            newEdge.set_vertices_pos(r1, r2)
            self.edges.append(newEdge)
            self.edges_vertices2index.update({newEdge.vertices:index})

            self.edges_num += 1

        print('{} edges added.'.format(self.edges_num))

    def build_faces(self):
        '''
        use ConvexHull to do triangulation of surface
        '''
        vertices = np.zeros( (self.vertices_num, 3) )
        for i in range(self.vertices_num):
            vertices[i][0] = self.vertices[i].x
            vertices[i][1] = self.vertices[i].y
            vertices[i][2] = self.vertices[i].z
        ConvexSurface = ConvexHull(vertices)
        surface_tri = ConvexSurface.simplices

        self.faces_num = surface_tri.shape[0]
        assert surface_tri.shape[1] == 3

        for i in range(self.faces_num):
            index = i
            v1 = surface_tri[i][0]
            v2 = surface_tri[i][1]
            v3 = surface_tri[i][2]

            r1 = vertices[v1]
            r2 = vertices[v2]
            r3 = vertices[v3]

            ort = int(np.sign( np.dot(r1, np.cross(r2, r3)) ))
            if ort == 1:
                newFace = Face(index, v1, v2, v3)
                newFace.set_vertices_pos(r1, r2, r3)
                self.faces.append(newFace)
                self.faces_vertices2index.update({newFace.vertices:index})
            else:
                newFace = Face(index, v3, v2, v1)
                newFace.set_vertices_pos(r3, r2, r1)
                self.faces.append(newFace)
                self.faces_vertices2index.update({newFace.vertices:index})

        print('{} faces added.'.format(self.faces_num))

    def build_dihedrals(self):
        '''
        build dihedrals according to faces
        '''
        self.dihedrals_num = 0
        dihedrals_set = set()

        for face in self.faces:
            # first common edge (v1 -> v2)
            edge1 = frozenset([face.v1, face.v2])

            adjacent_face = None
            for face_tmp in self.faces:
                if edge1.issubset(face_tmp.vertices) and (face_tmp.index != face.index):
                    adjacent_face = face_tmp
                    break
            if adjacent_face == None:
                raise Exception("Cannot find the adjacent triangle.")

            # v3 -> (v1 -> v2) -> v4
            v1 = face.v3
            v2 = face.v1
            v3 = face.v2
            v4 = next(iter(adjacent_face.vertices.difference(edge1)))
            dihedral_index = frozenset([v1,v2,v3,v4])
            if dihedral_index not in dihedrals_set:
                index = self.dihedrals_num

                r1 = self.vertices[v1].r
                r2 = self.vertices[v2].r
                r3 = self.vertices[v3].r
                r4 = self.vertices[v4].r

                newDihedral = Dihedral(index, v1, v2, v3, v4)
                newDihedral.set_vertices_pos(r1, r2, r3, r4)
                self.dihedrals.append(newDihedral)
                self.dihedrals_vertices2index.update({newDihedral.vertices:index})

                # edge opposite vertices info
                edge_index = self.edges_vertices2index[frozenset([v2,v3])]
                self.edges[edge_index].opposite_vertices.add(v1)
                self.edges[edge_index].opposite_vertices.add(v4)

                self.dihedrals_num += 1
                dihedrals_set.add(dihedral_index)

            # second common edge (v2 -> v3)
            edge2 = frozenset([face.v2, face.v3])

            adjacent_face = None
            for face_tmp in self.faces:
                if edge2.issubset(face_tmp.vertices) and (face_tmp.index != face.index):
                    adjacent_face = face_tmp
                    break
            if adjacent_face == None:
                raise Exception("Cannot find the adjacent triangle.")

            # v1 -> (v2 -> v3) -> v4
            v1 = face.v1
            v2 = face.v2
            v3 = face.v3
            v4 = next(iter(adjacent_face.vertices.difference(edge2)))
            dihedral_index = frozenset([v1,v2,v3,v4])
            if dihedral_index not in dihedrals_set:
                index = self.dihedrals_num

                r1 = self.vertices[v1].r
                r2 = self.vertices[v2].r
                r3 = self.vertices[v3].r
                r4 = self.vertices[v4].r

                newDihedral = Dihedral(index, v1, v2, v3, v4)
                newDihedral.set_vertices_pos(r1, r2, r3, r4)
                self.dihedrals.append(newDihedral)
                self.dihedrals_vertices2index.update({newDihedral.vertices:index})

                # edge opposite vertices info
                edge_index = self.edges_vertices2index[frozenset([v2,v3])]
                self.edges[edge_index].opposite_vertices.add(v1)
                self.edges[edge_index].opposite_vertices.add(v4)

                self.dihedrals_num += 1
                dihedrals_set.add(dihedral_index)

            # third common edge (v3 -> v1)
            edge3 = frozenset([face.v3, face.v1])

            adjacent_face = None
            for face_tmp in self.faces:
                if edge3.issubset(face_tmp.vertices) and (face_tmp.index != face.index):
                    adjacent_face = face_tmp
                    break
            if adjacent_face == None:
                raise Exception("Cannot find the adjacent triangle.")

            # v2 -> (v3 -> v1) -> v4
            v1 = face.v2
            v2 = face.v3
            v3 = face.v1
            v4 = next(iter(adjacent_face.vertices.difference(edge3)))
            dihedral_index = frozenset([v1,v2,v3,v4])
            if dihedral_index not in dihedrals_set:
                index = self.dihedrals_num

                r1 = self.vertices[v1].r
                r2 = self.vertices[v2].r
                r3 = self.vertices[v3].r
                r4 = self.vertices[v4].r

                newDihedral = Dihedral(index, v1, v2, v3, v4)
                newDihedral.set_vertices_pos(r1, r2, r3, r4)
                self.dihedrals.append(newDihedral)
                self.dihedrals_vertices2index.update({newDihedral.vertices:index})

                # edge opposite vertices info
                edge_index = self.edges_vertices2index[frozenset([v2,v3])]
                self.edges[edge_index].opposite_vertices.add(v1)
                self.edges[edge_index].opposite_vertices.add(v4)

                self.dihedrals_num += 1
                dihedrals_set.add(dihedral_index)

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
        self.r = np.array([x,y,z])

        # neighbors info
        self.neighbor_vertices = set()
        self.neighbor_edges = set()
        self.neighbor_faces = set()

        # properties
        self.type = 0
        self.normal = None
        self.mean_curvature = None
        self.gaussian_curvature = None
        self.curvature_1 = None
        self.curvature_2 = None
        self.area = None

    def set_pos(self, r):
        self.r = r

    def set_type(self, vtype):
        self.type = vtype

class Edge:
    '''
    represent an edge data (1D)
    '''
    def __init__(self, index, v1, v2):
        self.index = index
        self.v1 = v1
        self.v2 = v2
        self.vertices = frozenset([v1, v2])

        # neighbors info
        self.opposite_vertices = set()

        # properties
        self.type = 0

        self.r1 = None
        self.r2 = None

    def set_vertices_pos(self, r1, r2):
        self.r1 = r1
        self.r2 = r2

    def set_type(self, etype):
        self.type = etype

    def get_length(self):
        '''
        return edge length
        '''
        return np.linalg.norm(self.r1 - self.r2)

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

        self.r1 = None
        self.r2 = None
        self.r3 = None

        # properties
        self.type = 0
        self.area = None
        self.volume = None
        self.normal = None

    def set_vertices_pos(self, r1, r2, r3):
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3

    def set_type(self, ftype):
        self.type = ftype

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

        self.type = 0

        self.r1 = None
        self.r2 = None
        self.r3 = None
        self.r4 = None

    def set_vertices_pos(self, r1, r2, r3, r4):
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.r4 = r4

    def set_type(self, dtype):
        self.type = dtype

    def get_dihedral_angle(self):
        '''
        return dihedral angle

        following espresso convention

        v1 -> (v2 -> v3) -> v4
        '''
        r12 = self.r2 - self.r1
        r23 = self.r3 - self.r2
        r34 = self.r4 - self.r3

        n1 = np.cross(r12, r23)
        n2 = np.cross(r23, r34)

        n1_norm = np.linalg.norm(n1)
        n2_norm = np.linalg.norm(n2)

        n1_hat = n1 / n1_norm
        n2_hat = n2 / n2_norm

        cosphi = np.dot(n1_hat, n2_hat)

        phi = np.arccos(cosphi)

        return phi

class ElasticObjectType:
    '''
    Define properties of an elastic object
    '''
    def __init__(self, nodes_file="", triangles_file="", rescale=(1.0, 1.0, 1.0), ks_A=0.0, ks_B=0.0, kb_A=0.0, kb_B=0.0 ):
        # associated mesh data
        self.mesh = Mesh(nodes_file, triangles_file, rescale)

        # elastic parameters
        self.ks_A = ks_A
        self.ks_B = ks_B

        self.kb_A = kb_A
        self.kb_B = kb_B

        # bonded interactions
        self.stretching_interactions = []
        self.bending_interactions = []
        self.local_area_conservation = []
        self.global_area_conservation = []
        self.global_volume_conservation = []

    def initialize(self):
        '''
        initialization
        '''
        # initialize mesh data
        self.mesh.initialize()

        # initialize interactions
        self.define_bonded_interactions()

    def define_bonded_interactions(self):
        '''
        setup all bonded interactions
        '''
        # stretching
        self.define_stretching_interactions()

        # bending
        self.define_bending_interactions()

    def define_stretching_interactions(self):
        '''
        define stretching interactions for each edge
        '''
        # stretching
        for edge in self.mesh.edges:
            v1 = edge.v1
            v2 = edge.v2

            if edge.type == 0:
                # AA-bond
                ks_AA = self.ks_A / 2.0
                l0 = self.mesh.avg_edge_length
                hb_AA = espressomd.interactions.HarmonicBond(k=ks_AA, r_0=l0, r_cut=2*l0)
                self.stretching_interactions.append( (hb_AA, [v1,v2]) )
            elif edge.type == 1:
                # BB-bond
                ks_BB = self.ks_B / 2.0
                l0 = self.mesh.avg_edge_length
                hb_BB = espressomd.interactions.HarmonicBond(k=ks_BB, r_0=l0, r_cut=2*l0)
                self.stretching_interactions.append( (hb_BB, [v1,v2]) )
            elif edge.type == 2:
                # AB-bond
                ks_AB = self.ks_A*self.ks_B / (self.ks_A + self.ks_B)
                l0 = self.mesh.avg_edge_length
                hb_AB = espressomd.interactions.HarmonicBond(k=ks_AB, r_0=l0, r_cut=2*l0)
                self.stretching_interactions.append( (hb_AB, [v1,v2]) )
            else:
                raise Exception("Unrecognized edge type.")

    def define_bending_interactions(self):
        '''
        define bending interactions for each dihedral
        '''
        # bending
        for dihedral in self.mesh.dihedrals:
            v1 = dihedral.v1
            v2 = dihedral.v2
            v3 = dihedral.v3
            v4 = dihedral.v4

            if dihedral.type == 0:
                # AA-dihedral
                kb_AA = self.kb_A
                phi0 = np.pi
                #phi0 = dihedral.get_dihedral_angle()

                dihedral_AA = espressomd.interactions.Dihedral(bend=kb_AA, mult=1, phase=phi0)
                self.bending_interactions.append( (dihedral_AA, [v1,v2,v3,v4]) )
            elif dihedral.type == 1:
                # BB-dihedral
                kb_BB = self.kb_B
                phi0 = np.pi
                #phi0 = dihedral.get_dihedral_angle()

                dihedral_BB = espressomd.interactions.Dihedral(bend=kb_BB, mult=1, phase=phi0)
                self.bending_interactions.append( (dihedral_BB, [v1,v2,v3,v4]) )
            elif dihedral.type == 2:
                # AB-dihedral
                kb_AB = (self.kb_A + self.kb_B) / 2.0
                phi0 = np.pi
                #phi0 = dihedral.get_dihedral_angle()

                dihedral_AB = espressomd.interactions.Dihedral(bend=kb_AB, mult=1, phase=phi0)
                self.bending_interactions.append( (dihedral_AB, [v1,v2,v3,v4]) )
            else:
                raise Exception("Unrecognized dihedral type.")

    def define_local_area_conservation(self):
        '''
        define constrained forces which preserve the local area
        '''
        pass

    def define_global_area_conservation(self):
        '''
        define constrained forces which preserve the total area
        '''
        pass

    def define_global_volume_conservation(self):
        '''
        define constrained forces which preserve the total volume
        '''
        pass

class ElasticObject:
    '''
    represent an elastic object
    '''
    def __init__(self, system, elastic_object_type, object_id=0, particle_type_A=0, particle_type_B=1, translate=(0.0, 0.0, 0.0), rotate=(0.0, 0.0, 0.0)):
        # acctually create elastic object in espresso system
        self.system = system
        self.elastic_object_type = elastic_object_type
        self.mesh = elastic_object_type.mesh
        self.object_id = object_id
        self.particle_type_A = particle_type_A
        self.particle_type_B = particle_type_B
        self.index_offset = self.system.part.highest_particle_id + 1
    
    def initialize(self):
        '''
        initialization
        '''
        # create particles
        self.create_particles()

        # create bonded interactions
        #self.create_bonded_interactions()

    def create_particles(self):
        '''
        create espresso particles
        '''
        for vertex in self.mesh.vertices:
            particle_pos = vertex.r
            particle_type = self.particle_type_A if vertex.type == 0 else self.particle_type_B
            particle_id = self.index_offset + vertex.index
            self.system.part.add(pos=particle_pos, type=particle_type, id=particle_id)

    def create_bonded_interactions(self):
        '''
        create bonded interactions
        '''
        # stretching
        self.create_stretching_interactions()

        # bending
        self.create_bending_interactions()

    def create_stretching_interactions(self):
        '''
        create stretching interactions
        '''
        for entry in self.elastic_object_type.stretching_interactions:
            interaction = entry[0]
            partners = entry[1]
            assert len(partners) == 2

            p0 = partners[0]
            p1 = partners[1]

            self.system.bonded_inter.add(interaction)
            self.system.part[p0].add_bond( (interaction, p1) )
        
        print("stretching: {} bond interactions added.".format(len(self.elastic_object_type.stretching_interactions)))

    def create_bending_interactions(self):
        '''
        create bending interactions
        '''
        for entry in self.elastic_object_type.bending_interactions:
            interaction = entry[0]
            partners = entry[1]
            assert len(partners) == 4

            p0 = partners[0]
            p1 = partners[1]
            p2 = partners[2]
            p3 = partners[3]

            self.system.bonded_inter.add(interaction)
            self.system.part[p1].add_bond( (interaction, p0, p2, p3) )
        
        print("bending: {} dihedral interactions added.".format(len(self.elastic_object_type.bending_interactions)))

    def output_vtk_point_data_scalar(self, filename, dataname, datatype):
        '''
        write per-point scalar data
        '''
        with open(filename, 'a') as vtkfile:
            vtkfile.write("SCALARS {} {} 1\n".format(dataname, datatype))
            vtkfile.write("LOOKUP_TABLE default\n")
            for vertex in self.mesh.vertices:
                vtkfile.write('{}\n'.format(getattr(vertex, dataname)))

    def output_vtk_point_data_vector(self, filename, dataname, datatype):
        '''
        write per-point vector data
        '''
        with open(filename, 'a') as vtkfile:
            vtkfile.write("VECTORS {} {}\n".format(dataname, datatype))
            for vertex in self.mesh.vertices:
                vector = getattr(vertex, dataname)
                vtkfile.write('{} {} {}\n'.format(vector[0], vector[1], vector[2]))

    def output_vtk_cell_data_scalar(self, filename, dataname, datatype):
        '''
        write per-cell scalar data
        '''
        with open(filename, 'a') as vtkfile:
            vtkfile.write("SCALARS {} {} 1\n".format(dataname, datatype))
            vtkfile.write("LOOKUP_TABLE default\n")
            for face in self.mesh.faces:
                vtkfile.write('{}\n'.format(getattr(face, dataname)))

    def output_vtk_data(self, filename, output_attributes):
        '''
        write current configuration of elastic object
        '''

        with open(filename, 'w') as vtkfile:
            # Header
            vtkfile.write("# vtk DataFile Version 3.0\n")
            # Title
            vtkfile.write("Elastic Object {:d} Data\n".format(self.object_id))
            # Data type
            vtkfile.write("ASCII\n")
            # Geometry
            vtkfile.write("DATASET POLYDATA\n")
            # Points
            vtkfile.write("POINTS {:d} float\n".format(self.mesh.vertices_num))
            for vertex in self.mesh.vertices:
                # retrive current vertex position from espresso
                particle_id = vertex.index + self.index_offset
                particle_pos = self.system.part[particle_id].pos
                #particle_pos = self.system.part[particle_id].pos_folded
                vtkfile.write("{:f} {:f} {:f}\n".format(particle_pos[0], particle_pos[1], particle_pos[2]))
            # Vertices
            #vtkfile.write("VERTICES {:d} {:d}\n".format(self.mesh.vertices_num, self.mesh.vertices_num*2))
            #for vertex in self.mesh.vertices:
            #    # corresponding to point ids, no offset needed
            #    vtkfile.write("1 {:d}\n".format(vertex.index))
            # Lines
            #vtkfile.write("LINES {:d} {:d}\n".format(self.mesh.edges_num, self.mesh.edges_num*3))
            #for edge in self.mesh.edges:
            #    # corresponding to point ids, no offset needed
            #    v1 = edge.v1
            #    v2 = edge.v2
            #    vtkfile.write("2 {:d} {:d}\n".format(v1, v2))
            # Triangles
            vtkfile.write("POLYGONS {:d} {:d}\n".format(self.mesh.faces_num, self.mesh.faces_num*4))
            for face in self.mesh.faces:
                # correspond to point ids, no offset needed
                v1 = face.v1
                v2 = face.v2
                v3 = face.v3
                vtkfile.write("3 {:d} {:d} {:d}\n".format(v1, v2, v3))
            
        # additional attributes

        # point data - scalar (default)
        with open(filename, 'a') as vtkfile:
            vtkfile.write("POINT_DATA {}\n".format(self.mesh.vertices_num))
            vtkfile.write("SCALARS {} {} 1\n".format("index", "int"))
            vtkfile.write("LOOKUP_TABLE default\n")
            for vertex in self.mesh.vertices:
                vtkfile.write("{:d}\n".format(vertex.index))

        # vertex type
        if 'type' in output_attributes:
            self.output_vtk_point_data_scalar(filename, 'type', 'int')
        
        # vertex mean curvature
        if 'mean_curvature' in output_attributes:
            self.output_vtk_point_data_scalar(filename, 'mean_curvature', 'float')

        # vertex gaussian curvature
        if 'gaussian_curvature' in output_attributes:
            self.output_vtk_point_data_scalar(filename, 'gaussian_curvature', 'float')
        
        # vertex principal curvatures
        if 'principal_curvatures' in output_attributes:
            self.output_vtk_point_data_scalar(filename, 'curvature_1', 'float')
            self.output_vtk_point_data_scalar(filename, 'curvature_2', 'float')
        
        # vertex area
        if 'area' in output_attributes:
            self.output_vtk_point_data_scalar(filename, 'area', 'float')

        # vertex normal
        if 'normal' in output_attributes:
            self.output_vtk_point_data_vector(filename, 'normal', 'float')
