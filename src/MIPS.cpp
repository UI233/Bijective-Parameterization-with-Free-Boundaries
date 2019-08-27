#include "MIPS.h"
#include "Eigen/SparseCholesky"
#include "Eigen/Core"
#include <iostream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <limits>
#include <queue>
#include "OpenMesh/Core/IO/MeshIO.hh"

#define TEST
// the hashmap for neighbor search
class myHash {
public:
    static constexpr size_t maxsize = 19260817u;
    size_t operator()(const OpenMesh::Vec2i& idx) const {
        return maxsize * idx[0] + idx[1];
    }
};

static std::unordered_multimap<OpenMesh::Vec2i, int, myHash> neighbors;

template<typename T>
auto sq(const T& v) -> decltype(v * v) {
    return v * v;
}

// Project the point onto the unit square from a circle
OpenMesh::Vec2d project2Circle(double theta) {
    return OpenMesh::Vec2d(cos(theta), sin(theta));
}

// Get the initial parameterization of the mesh via
// minimizing the energy E = 1/2 \sum{|| p_i - p_j ||} for the inner point while the boundary points are projected onto the unit circle (doubleer Method)
// It is problemetic I think.
std::vector<OpenMesh::Vec2f> init(Mesh& mesh) {
    std::set<Mesh::VertexHandle> boundaries;
    std::vector<OpenMesh::Vec2f> res;
    std::vector<int> inner;
    size_t vertices_sz = 0;
    size_t inner_count = -1;
    constexpr int coord = 2;

    // select the boundary
    for (auto cit = mesh.vertices_begin(); cit != mesh.vertices_end(); ++cit) {
        if (mesh.is_boundary(*cit)) {
            boundaries.insert(*cit);
        }
        else inner_count++;
        inner.push_back(inner_count);
        ++vertices_sz;
    }
    res.resize(vertices_sz);
    Eigen::SparseMatrix<double> sparse(coord * (vertices_sz - boundaries.size()), coord * (vertices_sz - boundaries.size()));

    auto getNext = [&mesh] (OpenMesh::VertexHandle st) {
        for (auto f: mesh.vf_range(st)) {
            OpenMesh::VertexHandle vs[3];
            int i = 0, iv = 0;
            for (auto e_itr = mesh.fe_ccwbegin(f); e_itr.is_valid(); ++e_itr) {
                auto e = *e_itr;
                if (mesh.is_boundary(e)) {
                    auto s = mesh.from_vertex_handle(mesh.halfedge_handle(e, 0));
                    auto t = mesh.to_vertex_handle(mesh.halfedge_handle(e, 0));
                    
                    if (s == st)
                        return t;
                }
            }
        }

        return OpenMesh::VertexHandle(-1);
    };
    // compute the boundary points
    if (boundaries.size() > 0) {
        std::vector<bool> visited;
        visited.resize(vertices_sz);
        std::fill(visited.begin(), visited.end(), 0);
        double length = 0.0;

        // compute the total length
        for (auto st = *boundaries.begin(); !visited[st.idx()];) {
            visited[st.idx()] = true;
            auto nx = getNext(st);
            length += (mesh.point(nx) - mesh.point(st)).norm();
            st = nx;
        }

        // compute the boundary point
        std::fill(visited.begin(), visited.end(), 0);
        double angle = 0.0;
        int cnt = 0;
        auto bst = *boundaries.begin();
        for (auto st = *boundaries.begin(); !visited[st.idx()];) {
            visited[st.idx()] = true;
            ++cnt;
            res[st.idx()] = project2Circle(angle);
            auto nx = getNext(st);
            float len = (mesh.point(nx) - mesh.point(st)).norm();
            angle += 2.0 * M_PI * len / length;
            st = nx;
        }
    }
    
    if (boundaries.size() != vertices_sz) {
        // initialize the boundary point and the column
        Eigen::VectorXd B(coord * (vertices_sz - boundaries.size())); 
        B.setZero();
        for (auto bd: boundaries) {
            auto point = mesh.point(bd);
            auto p = res[bd.idx()];
            res[bd.idx()] = p;
            for (auto it = mesh.vv_begin(bd); it.is_valid(); ++it) 
                if (!mesh.is_boundary(*it)) {
                    int ii = inner[(*it).idx()];
                    auto dis = (point  - mesh.point(*it)).norm();
                    sparse.coeffRef(coord * ii,  coord * ii) += 1.0 / 1.0; // / dis;
                    sparse.coeffRef(coord * ii + 1,  coord * ii + 1) += 1.0 / 1.0; // / dis;
                    B[coord * ii] += p[0] / 1.0; // / dis;
                    B[coord * ii + 1] += p[1] / 1.0; // / dis;
                }
        }
        // build the sparse matrix
        for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i) {
            if(!boundaries.count(*i)) {
                for (auto j = mesh.vv_iter(*i); j.is_valid(); ++j) {
                    if (!boundaries.count(*j) && *i != *j) {
                        int ii = inner[(*i).idx()], ij = inner[(*j).idx()];
                        double dis = (mesh.point(*i) - mesh.point(*j)).norm();
                        for (int cnt = 0; cnt < coord; ++cnt) {
                            sparse.coeffRef(coord * ii + cnt, coord * ij + cnt) -= 0.5 / 1.0; // / dis;
                            sparse.coeffRef(coord * ij + cnt, coord * ii + cnt) -= 0.5 / 1.0; // / dis;
                            sparse.coeffRef(coord * ii + cnt, coord * ii + cnt) += 1.0 / 1.0; // / dis;
                        }
                    }
                }
            }
        }

        // solve the system to get the initial points set
        Eigen::SimplicialLDLT<decltype(sparse)> solver;    
        auto X = solver.compute(sparse).solve(B);
        int i = 0;
        for (int j = 0; j < vertices_sz; ++j)
            if (!boundaries.count(Mesh::VertexHandle(j))){
                res[j][0] = X[i]; 
                res[j][1] = X[i + 1]; 
                i += 2;
            }
    }
    return std::move(res);
}


// initialize the hash map for neighbor search
void initHashMap(const Mesh& mesh, const Eigen::VectorXd& x, double epsilon) {
    neighbors.clear();
    for (int i = 0; i < x.size(); i += 2) {
        OpenMesh::VertexHandle handle(i / 2);
        if (mesh.is_boundary(handle)) {
            OpenMesh::Vec2i idx(x[i] / epsilon, x[i + 1] / epsilon);
            neighbors.insert(std::make_pair(idx, i / 2));
        }
    }
}

double getDistance(const Mesh& mesh, const Eigen::VectorXd& params, const OpenMesh::EdgeHandle& e, const OpenMesh::VertexHandle& p) {
    int st = mesh.to_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
    int end = mesh.from_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
    OpenMesh::Vec2d pst(params[2 * st], params[2 * st + 1]); 
    OpenMesh::Vec2d pend(params[2 * end], params[2 * end + 1]); 
    OpenMesh::Vec2d edge = pend - pst;

    OpenMesh::Vec2d point(params[p.idx() * 2], params[p.idx() * 2 + 1]);
    double projected = OpenMesh::dot(point - pst, edge) / edge.norm();
    double dis = sqrtf(OpenMesh::dot(point - pst, point - pst) - projected);

    if (OpenMesh::dot(point - pend, pst - pend) < 0.0) 
        dis = (point - pend).norm();
    else if (OpenMesh::dot(point - pst, pend - pst) < 0.0) 
        dis = (point - pst).norm();

    return dis;
}

// Get the energy E_D + E_B for current parameterization and ε
double getEnergy(const Mesh& mesh, const Eigen::VectorXd& params, double epsilon, bool display = false) {
    // update the hashmap for computing the energy and gradient
    initHashMap(mesh, params, epsilon);
    double energy = 0.0;
    using OpenMesh::dot;
    // distortion term
    for (auto f: mesh.all_faces()) {
        int i = 0;
        Mesh::VertexHandle vs[3];
        for (auto v: mesh.fv_range(f))
            vs[i++] = v;
        OpenMesh::Vec3f u0(params[2 * vs[0].idx()], params[2 * vs[0].idx() + 1], 0.0);
        OpenMesh::Vec3f ue1 = OpenMesh::Vec3f(params[2 * vs[1].idx()], params[2 * vs[1].idx() + 1], 0.0) - u0;
        OpenMesh::Vec3f ue2 = OpenMesh::Vec3f(params[2 * vs[2].idx()], params[2 * vs[2].idx() + 1], 0.0) - u0;

        OpenMesh::Vec3f p0 = mesh.point(vs[0]);
        OpenMesh::Vec3f pe1 = mesh.point(vs[1]) - p0;
        OpenMesh::Vec3f pe2 = mesh.point(vs[2]) - p0;

        // calculate the area of corresponding triangle
        double su = 0.5f * OpenMesh::cross(ue1, ue2).norm();
        double sp = 0.5f * OpenMesh::cross(pe1, pe2).norm();

        energy += (1 + (sp * sp) / (su * su)) * (sq(ue2.norm()) * sq(pe1.norm()) + sq(ue1.norm()) * sq(pe2.norm()) - 2.0f * dot(ue2, ue1) * dot(pe2, pe1)) / (4.0f * sp);
    }

    double tseng = 0.0;
    // the barrier term
    for (auto e: mesh.all_edges())
        if (mesh.is_boundary(e)) {
            int st = mesh.to_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
            int end = mesh.from_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
            // OpenMesh::Vec2d pst(params[2 * st], params[2 * st + 1]); 
            // OpenMesh::Vec2d pend(params[2 * end], params[2 * end + 1]); 
            // OpenMesh::Vec2d edge = pend - pst;

            // The range for neighborhood searching
            int x_st = std::min(params[2 * st], params[2 * end]) / epsilon - 1;
            int x_end = std::max(params[2 * st], params[2 * end]) / epsilon + 1;
            int y_st = std::min(params[2 * st + 1], params[2 * end + 1]) / epsilon - 1;
            int y_end = std::max(params[2 * st + 1], params[2 * end + 1]) / epsilon + 1;
            // conduct neighbor search for neighbor boundaries points
            for (int x = x_st; x <= x_end; ++x)
                for (int y = y_st; y <= y_end; ++y) {
                    auto itr_pr = neighbors.equal_range(OpenMesh::Vec2i(x, y));
                    for (auto itr = itr_pr.first; itr != itr_pr.second; ++itr)
                        if (itr->second != st && itr->second != end) {
                            // OpenMesh::Vec2d point(params[itr->second * 2], params[itr->second * 2 + 1]);
                            // double projected = OpenMesh::dot(point - pst, edge) / edge.norm();
                            // double dis = sqrtf(OpenMesh::dot(point - pst, point - pst) - projected * projected);
                            // if (OpenMesh::dot(point - pend, pst - pend) < 0.0) 
                            //     dis = (point - pend).norm();
                            // else if (OpenMesh::dot(point - pst, pend - pst) < 0.0) 
                            //     dis = (point - pst).norm();
                            // std::cout << "Dis: " << dis << std::endl;
                            double dis = getDistance(mesh, params, e, OpenMesh::VertexHandle(itr->second));
                            if (dis == 0.0)
                                return std::numeric_limits<double>::infinity();
                            double seng = std::max(0.0, epsilon / dis - 1.0);
                            tseng += seng * seng;
                        }
                }
        }

    energy += tseng;
    
    // double test_seng = 0.0;
    // for (auto e: mesh.all_edges()) 
    //     if (mesh.is_boundary(e)) {
    //         int st = mesh.to_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
    //         int end = mesh.from_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
    //         int cnt = 0;
    //         for (auto v: mesh.all_vertices())
    //             if (v.idx() != st  && v.idx() != end && mesh.is_boundary(v)) {
    //                 ++cnt;
    //                 OpenMesh::Vec2d pst(params[2 * st], params[2 * st + 1]); 
    //                 OpenMesh::Vec2d pend(params[2 * end], params[2 * end + 1]); 
    //                 OpenMesh::Vec2d edge = pend - pst;
    //                 OpenMesh::Vec2d point(params[v.idx() * 2], params[v.idx() * 2 + 1]);
    //                 double projected = OpenMesh::dot(point - pst, edge) / edge.norm();
    //                 double dis = sqrtf(OpenMesh::dot(point - pst, point - pst) - projected * projected);
    //                 if (OpenMesh::dot(point - pend, pst - pend) < 0.0) 
    //                     dis = (point - pend).norm();
    //                 else if (OpenMesh::dot(point - pst, pend - pst) < 0.0) 
    //                     dis = (point - pst).norm();

    //                 double seng = std::max(0.0, epsilon / dis - 1.0);
    //                 test_seng += seng * seng;
    //             }
    //     }
    
    if(display)
        std::cout << "Rport:\nBoundary Barrier: [" << tseng << "]" << std::endl <<
            "Distortion: [" << energy - tseng << "]\n"  << std::endl << 
            "Total: " << energy << std::endl;
    return energy;
}

// calculate the gradient of the energy function
auto calcGrad(const Mesh& mesh, const Eigen::VectorXd& x, double epsilon) {
    // update the hashmap for computing the energy and gradient
    initHashMap(mesh, x, epsilon);
    // the energy term
    Eigen::VectorXd gradx(x.size());
    gradx.setZero();
    // the gradient for distortion term
    for (auto f: mesh.all_faces()) {
        Mesh::VertexHandle vs[3];
        int i = 0;
        for (auto itr = mesh.cfv_ccwbegin(f); itr != mesh.cfv_ccwend(f); itr++)
            vs[i++] = *itr;
        OpenMesh::Vec3d u0(x[2 * vs[0].idx()], x[2 * vs[0].idx() + 1], 0.0);
        OpenMesh::Vec3d u1(x[2 * vs[1].idx()], x[2 * vs[1].idx() + 1], 0.0);
        OpenMesh::Vec3d u2(x[2 * vs[2].idx()], x[2 * vs[2].idx() + 1], 0.0);
        OpenMesh::Vec3d ue1 = u1 - u0;
        OpenMesh::Vec3d ue2 = u2 - u0;

        OpenMesh::Vec3d p0(mesh.point(vs[0])[0], mesh.point(vs[0])[1], mesh.point(vs[0])[2]);
        OpenMesh::Vec3d p1(mesh.point(vs[1])[0], mesh.point(vs[1])[1], mesh.point(vs[1])[2]);
        OpenMesh::Vec3d p2(mesh.point(vs[2])[0], mesh.point(vs[2])[1], mesh.point(vs[2])[2]);
        OpenMesh::Vec3d pe1 = p1 - p0;
        OpenMesh::Vec3d pe2 = p2 - p0;

        // calculate the area of corresponding triangle
        double su = 0.5 * OpenMesh::cross(ue1, ue2).norm();
        if (su == 0.0) {
            gradx[0] = std::numeric_limits<double>::infinity();
            return std::move(gradx);
        }
        double sp = 0.5 * OpenMesh::cross(pe1, pe2).norm();

        double cot0 = OpenMesh::dot(p1 - p0, p2 - p0) / OpenMesh::cross(p1 - p0, p2 - p0).norm();
        double cot1 = OpenMesh::dot(p2 - p1, p0 - p1) / OpenMesh::cross(p2 - p1, p0 - p1).norm();
        double cot2 = OpenMesh::dot(p1 - p2, p0 - p2) / OpenMesh::cross(p1 - p2, p0 - p2).norm();
        
        double term1 = (1 + (sp * sp) / (su * su));
        double term2 = (sq(ue2.norm()) * sq(pe1.norm()) + sq(ue1.norm()) * sq(pe2.norm()) - 2.0f * dot(ue2, ue1) * dot(pe2, pe1)) / (4.0f * sp);

        {
            // dE / du0
            OpenMesh::Vec2d height(-(u2 - u1)[1], (u2 - u1)[0]);
            if (OpenMesh::dot(OpenMesh::Vec3d(height[0], height[1], 0.0f), u0 - u1) < 0.0f)  {
                height = -height;
            }
            OpenMesh::Vec2d dterm1_du0 = - (sp* sp) / (su * su * su) * height;
            OpenMesh::Vec3d dterm2_du0 = -cot2 * u1 - cot1 * u2 + (cot1 + cot2) * u0;
            OpenMesh::Vec2d dE_du0 = dterm1_du0 * term2 + OpenMesh::Vec2d(dterm2_du0[0], dterm2_du0[1]) * term1; 
            gradx[vs[0].idx() * 2] += dE_du0[0];
            gradx[vs[0].idx() * 2 + 1] += dE_du0[1];
        }

        {
            // dE / du1
            OpenMesh::Vec2d height(-(u0 - u2)[1], (u0 - u2)[0]);
            if (OpenMesh::dot(OpenMesh::Vec3d(height[0], height[1], 0.0f), u1 - u0) < 0.0f)  {
                height = -height;
            }
            OpenMesh::Vec2d dterm1_du1 = - (sp* sp) / (su * su * su) * height;
            OpenMesh::Vec3d dterm2_du1 = -cot2 * u0 - cot0 * u2 + (cot0 + cot2) * u1;
            OpenMesh::Vec2d dE_du1 = dterm1_du1 * term2 + OpenMesh::Vec2d(dterm2_du1[0], dterm2_du1[1]) * term1; 
            gradx[vs[1].idx() * 2] += dE_du1[0];
            gradx[vs[1].idx() * 2 + 1] += dE_du1[1];
        }

        {
            // dE / du2
            OpenMesh::Vec2d height(-(u1 - u0)[1], (u1 - u0)[0]);
            if (OpenMesh::dot(OpenMesh::Vec3d(height[0], height[1], 0.0f), u2 - u0) < 0.0f)  {
                height = -height;
            }
            OpenMesh::Vec2d dterm1_du2 = -(sp* sp) / (su * su * su) * height;
            OpenMesh::Vec3d dterm2_du2 = -cot0 * u1 - cot1 * u0 + (cot0 + cot1) * u2;
            OpenMesh::Vec2d dE_du2 = dterm1_du2 * term2 + OpenMesh::Vec2d(dterm2_du2[0], dterm2_du2[1]) * term1; 
            gradx[vs[2].idx() * 2] += dE_du2[0];
            gradx[vs[2].idx() * 2 + 1] += dE_du2[1];
        }
    }

    int cnt = 0;
    // the singularity term
    for (auto e: mesh.all_edges())
        if (mesh.is_boundary(e)) {
            double seng = 0.0;
            int ui_idx = -1;

            int st = mesh.to_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
            int end = mesh.from_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
            OpenMesh::Vec2d pst(x[2 * st], x[2 * st + 1]);
            OpenMesh::Vec2d pend(x[2 * end], x[2 * end + 1]);
            OpenMesh::Vec2d edge = pend - pst;

            // The range for neighborhood searching
            int x_st = std::min(x[2 * st], x[2 * end]) / epsilon - 1;
            int x_end = std::max(x[2 * st], x[2 * end]) / epsilon + 1;
            int y_st = std::min(x[2 * st + 1], x[2 * end + 1]) / epsilon - 1;
            int y_end = std::max(x[2 * st + 1], x[2 * end + 1]) / epsilon + 1;
            // conduct neighbor search for neighbor boundaries points
            for (int xs = x_st; xs <= x_end; ++xs)
                for (int y = y_st; y <= y_end; ++y) {
                    auto itr_pr = neighbors.equal_range(OpenMesh::Vec2i(xs, y));
                    for (auto itr = itr_pr.first; itr != itr_pr.second; ++itr)
                        if (itr->second != st && itr->second != end) {
                            OpenMesh::Vec2d point(x[itr->second * 2], x[itr->second * 2 + 1]);
                            double projected = OpenMesh::dot(point - pst, edge) / edge.norm();
                            double dis = getDistance(mesh, x, e, OpenMesh::VertexHandle(itr->second));

                            if (dis == 0.0) {
                                gradx[0] = std::numeric_limits<double>::infinity();
                                return std::move(gradx);
                            }
                            seng = std::max(0.0, epsilon / dis - 1.0);
                            // calculate the gradient
                            ui_idx = itr->second;
                            if (ui_idx != -1 && seng != 0.0) {
                                OpenMesh::Vec2d grad_maxui;
                                OpenMesh::Vec2d grad_maxu1;
                                OpenMesh::Vec2d grad_maxu2;
                                OpenMesh::Vec2d ui(x[ui_idx * 2], x[ui_idx * 2 + 1]);
                                auto u1 = pst;
                                auto u2 = pend;
                                double dis = sqrt(OpenMesh::dot(ui - u1, ui - u1) - projected * projected);
                                
                                if (OpenMesh::dot(point - pend, pst - pend) < 0.0) {
                                    OpenMesh::Vec2d dir_ui = (point - pend).normalize();
                                    OpenMesh::Vec2d dir_u2 = -(point - pend).normalize();
                                    
                                    // dE / du2
                                    gradx[end * 2] += -2.0 * seng * epsilon * dir_u2[0] / (dis * dis);
                                    gradx[end * 2 + 1] +=  -2.0 * seng * epsilon * dir_u2[1] / (dis * dis);

                                    // dE / dui
                                    gradx[ui_idx * 2] += -2.0 * seng * epsilon * dir_ui[0] / (dis * dis);
                                    gradx[ui_idx * 2 + 1] += -2.0 * seng * epsilon * dir_ui[1] / (dis * dis);
                                }
                                else if (OpenMesh::dot(point - pst, pend - pst) < 0.0) {
                                    OpenMesh::Vec2d dir_ui = (point - pst).normalize();
                                    OpenMesh::Vec2d dir_u1 = -(point - pst).normalize();

                                    // dE / du1
                                    gradx[st * 2] += -2.0 * seng * epsilon * dir_u1[0] / (dis * dis);
                                    gradx[st * 2 + 1] +=  -2.0 * seng * epsilon * dir_u1[1] / (dis * dis);

                                    // dE / dui
                                    gradx[ui_idx * 2] += -2.0 * seng * epsilon * dir_ui[0] / (dis * dis);
                                    gradx[ui_idx * 2 + 1] += -2.0 * seng * epsilon * dir_ui[1] / (dis * dis);

                                }
                                else {
                                    double projected = OpenMesh::dot(ui - u1, u2 - u1) / edge.norm();
                                    // the unit vector othorgonal to the edge
                                    auto normal = ui - u1 - (u2 - u1).normalize() * projected;
                                    normal = normal.normalize();
                                    // dE / du2 u2 is the end point of this edge
                                    grad_maxu2 = -2 * seng * epsilon * projected / OpenMesh::dot(u2 - u1, u2 - u1) * normal / (dis * dis);
                                    gradx[end * 2] += grad_maxu2[0];
                                    gradx[end * 2 + 1] += grad_maxu2[1];
                                    // dE / du1 u1 is the start point of this edge
                                    projected = OpenMesh::dot(ui - u2, u1 - u2) / edge.norm();
                                    grad_maxu1 = -2 * seng * epsilon * projected / OpenMesh::dot(u2 - u1, u2 - u1) * normal / (dis * dis);
                                    gradx[st * 2] += grad_maxu1[0];
                                    gradx[st * 2 + 1] += grad_maxu1[1];
                                    // dE / dui ui is the point with maximal distance
                                    grad_maxui = -2 * seng * epsilon * normal / (dis * dis);
                                    gradx[ui_idx * 2] += grad_maxui[0];
                                    gradx[ui_idx * 2 + 1] += grad_maxui[1];
                                }
                            }
                        }
                }
        }
    return std::move(gradx);
}

// compute the average length of edges
double getAveLen(const Mesh& mesh, const Eigen::VectorXd& x) {
    int count = 0;
    double len = 0.0;
    for (auto e: mesh.all_edges()) {
        ++count;
        int st = mesh.to_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
        int end = mesh.from_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
        OpenMesh::Vec2d ev(x[2 * st] - x[2 *end], x[2 * st + 1] - x[2 * end + 1]);
        len += ev.norm();
    }

    return len / count; 
}



bool checkCreterion(const Mesh& mesh, const Eigen::VectorXd& x, double epsilon) {
    initHashMap(mesh, x, epsilon);
    for (auto f: mesh.all_faces()) {
        // get the vertices associated with this triangle
        OpenMesh::VertexHandle vs[3];
        OpenMesh::Vec2d param[3];
        OpenMesh::Vec3d pos3d[3];
        OpenMesh::Vec2d u1 = param[1] - param[0];
        OpenMesh::Vec2d u2 = param[2] - param[0];

        if (u2[1] * u1[0] == u2[0] * u1[1])
            return false;

        int count = 0;
        for (auto v: mesh.fv_range(f)) {
            vs[count] = v;
            pos3d[count] = mesh.point(v);
            param[count] = OpenMesh::Vec2d(x[v.idx() * 2], x[v.idx() * 2 + 1]);
            ++count;
        }
        // check non-flip and bijective
        OpenMesh::Vec2d e1 = (param[1] - param[0]) / (pos3d[1] - pos3d[0]).norm();
        double xc = OpenMesh::dot(pos3d[2] - pos3d[0], pos3d[1] - pos3d[0]) / (pos3d[1] - pos3d[0]).norm();
        double yc = sqrt(OpenMesh::dot(pos3d[2] - pos3d[0], pos3d[2] - pos3d[0]) - xc * xc);
        OpenMesh::Vec2d e2 =  (param[2] - param[0] - e1 * xc) / yc; 
        double x00 = e1[0];
        double x10 = e1[1];
        double x01 = e2[0];
        double x11 = e2[1];
        double det = x00 * x11 - x01 * x10;
        if (det < 0.0)
            return false;
    }

    // check the degeneration of boundary points
    for (auto e: mesh.all_edges())
        if (mesh.is_boundary(e)) {
            int st = mesh.from_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
            int end = mesh.to_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
            // The range for neighborhood searching
            int x_st = std::min(x[2 * st], x[2 * end]) / epsilon - 1;
            int x_end = std::max(x[2 * st], x[2 * end]) / epsilon + 1;
            int y_st = std::min(x[2 * st + 1], x[2 * end + 1]) / epsilon - 1;
            int y_end = std::max(x[2 * st + 1], x[2 * end + 1]) / epsilon + 1;
            // conduct neighbor search for neighbor boundaries points
            for (int xs = x_st; xs <= x_end; ++xs)
                for (int y = y_st; y <= y_end; ++y) {
                    auto itr_pr = neighbors.equal_range(OpenMesh::Vec2i(xs, y));
                    for (auto itr = itr_pr.first; itr != itr_pr.second; ++itr)
                        if (itr->second != st && itr->second != end) {
                            double dis = getDistance(mesh, x, e, OpenMesh::VertexHandle(itr->second));
                            if (dis == 0.0)
                                return false;
                        }
                    }
        }
    return true;
}

// get the step length for line search
double getStepLength(const Mesh& mesh, const Eigen::VectorXd& params, const Eigen::VectorXd& p, const Eigen::VectorXd& grad, const double& epsilon, double alpha0) {
    constexpr double c1 = 1e-4;
    double ene0 = getEnergy(mesh, params, epsilon);
    double enei = getEnergy(mesh, params + alpha0 * p, epsilon);
    double dphi = p.dot(grad);
    while (enei > ene0 + c1 * alpha0 * dphi || !checkCreterion(mesh, params + alpha0 * p, epsilon)) {
        alpha0 *= 0.50;
        enei = getEnergy(mesh, params + alpha0 * p, epsilon);
    }

    std::cout << "Step len: " << alpha0 << " dphi: " << dphi << std::endl;
    return alpha0;
    // if (enei < ene0 + c1 * alpha0 * dphi)
    //     return  alpha0;
    // double alpha = - dphi * alpha0 * alpha0 / (2.0f * (enei - ene0 - dphi * alpha0));
    // double enei_0 = enei; // energy_{i-1}
    // enei = getEnergy(mesh, params + alpha * p, epsilon);
    // while (!(enei < ene0 + c1 * alpha * dphi)) {
    //     // calculate the coeff of the cubic interpolation
    //     double dim0 = enei - ene0 - dphi * alpha;
    //     double dim1 = enei_0 - ene0 - dphi * alpha0;
    //     double numerator = 1.0 / (alpha0 * alpha0 * alpha * alpha * (alpha - alpha0));

    //     double a = (alpha0 * alpha0 * dim0 - alpha * alpha * dim1) * numerator;
    //     double b = (-alpha0 * alpha0 * alpha0 * dim0 + alpha * alpha * alpha * dim1) * numerator;
    //     // calculate the aplha
    //     alpha0 = alpha;
    //     alpha = (-b + sqrt(b * b - 3 * a *dphi)) / (3 * a);

    //     if (fabs(alpha - alpha0) / fabs(alpha0) < 0.11)
    //         alpha = alpha0 / 2.0f;
    //     // update 
    //     enei_0 = enei;
    //     enei = getEnergy(mesh, params + alpha * p, epsilon);
    // }
    // return alpha;
}

// To be tested
auto LBFGS(Mesh& mesh, Eigen::VectorXd x, const double& error) {
    static constexpr int sz = 7;
    auto test = [&error](const Eigen::VectorXd& x) -> bool {
        return x.norm() > error;
    };

    Eigen::VectorXd ss[sz], ys[sz];
    double rhos[sz], alphas[sz];
    double epsilon = getAveLen(mesh, x) * 0.25f;
    Eigen::VectorXd gradx = calcGrad(mesh, x, epsilon);
    int k = 0;
    int step = 0;
    do{
        #ifdef TEST
            Mesh omesh;
            for (int i = 0; i < x.size(); i += 2) {
                OpenMesh::Vec3f vert(x[i], x[i + 1], 0.0f);
                omesh.add_vertex(vert);
            }
            for (auto face: mesh.all_faces()) {
                int cnt = 0;
                OpenMesh::VertexHandle handles[10];
                for (auto v: mesh.fv_range(face)) 
                    handles[cnt++] = v;
                omesh.add_face(handles, cnt);
            }
            if(k % 10 == 0 || k == 1)
                OpenMesh::IO::write_mesh(omesh, std::to_string(k) + "test" + "_param" +  ".obj");
        #endif
        auto q = gradx;
        auto gradx0 = gradx;
        Eigen::VectorXd r(x.size());
        for (int i = k - 1; i >= std::max(0, k - sz); --i) {
            alphas[i % sz] = q.dot(ss[i % sz]) * rhos[i % sz];
            q -= alphas[i % sz] * ys[i % sz]; 
        }
        if (k == 0)
            r = q;
        else r = ss[(k - 1) % sz].dot(ys[(k - 1) % sz]) / ys[(k - 1) % sz].dot(ys[(k - 1) % sz]) * q;
        for (int i = std::max(0, k - sz); i < k; ++i) {
            auto beta = rhos[i % sz] * r.dot(ys[i %sz]);
            r += ss[i % sz] * (alphas[i % sz] - beta);
        }
        double step_len = getStepLength(mesh, x, -r, gradx0, epsilon, k == 0 ? 1e-5 : 1.0);
        auto energy = getEnergy(mesh, x, epsilon, true);
        std::cout << "Step: " << step++ << " Norm of Gradient: "<< gradx.norm() << std::endl;
        // update array for y, s, rho
        ss[k % sz] = step_len * -r;
        x += ss[k % sz];
        gradx = calcGrad(mesh, x, epsilon);
        ys[k % sz] = gradx - gradx0;
        rhos[k % sz] = 1.0 / ss[k % sz].dot(ys[k % sz]);
       
        ++k;
    }while(test(gradx)); // todo!

    return std::move(x);
}

auto gradientDescent(Mesh& mesh, Eigen::VectorXd x, const double& error) {
    int step = 0;
    double epsilon = 0.25f * getAveLen(mesh, x);
    do{
        auto energy = getEnergy(mesh, x, epsilon, true);
        auto gradx = calcGrad(mesh, x, epsilon); 
        std::cout << "Step: " << step++ << " Norm of Gradient: "<< gradx.norm() << std::endl;
        auto len = getStepLength(mesh, x, -gradx, gradx, epsilon, 1e-5);
        x -= len * gradx;
    }while(true);
    
    return std::move(x);
}

std::vector<OpenMesh::Vec2f> paratimization(Mesh& mesh, double error) {
    auto init_param = init(mesh); // the initial parameterization for the mesh
    Eigen::VectorXd x(init_param.size() * 2);
    for (int i = 0; i < init_param.size(); ++i) {
        x[2 * i] = init_param[i][0];
        x[2 * i + 1] = init_param[i][1];
    }

    auto epsilon = getAveLen(mesh, x) * 0.25;
    std::cerr << checkCreterion(mesh, x, epsilon) << std::endl;
    x = LBFGS(mesh, x, 0.0005); // optimize the parameterization using Quasi Newton Method
    // x = gradientDescent(mesh, x, 0.05);

    for (int i = 0; i < init_param.size(); ++i) {
        init_param[i][0] = x[2 * i];
        init_param[i][1] = x[2 * i + 1];
    }

    return std::move(init_param);
}

void test(Mesh& mesh, const std::vector<OpenMesh::Vec2f>& params) {
    Eigen::VectorXd x(params.size() * 2);
    for (int i = 0; i < params.size(); ++i) {
        x[2 * i] = params[i][0]; 
        x[2 * i + 1] = params[i][1]; 
    }
    auto epsilon = getAveLen(mesh, x) * 0.25f;
    auto gradx = calcGrad(mesh, x, epsilon);
    getEnergy(mesh, x, epsilon, true);
    std::cout << "Log: [Gradient]\n";
    for (int i = 0; i < params.size(); ++i) {
        std::cout << "(" << gradx[i * 2] << ", " << gradx[i * 2 + 1] << ")" << std::endl;
    }
}