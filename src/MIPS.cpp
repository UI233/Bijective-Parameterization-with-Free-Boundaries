#include "MIPS.h"
#include "Eigen/SparseCholesky"
#include "Eigen/Core"
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <limits>
#include <queue>

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
OpenMesh::Vec2f project2Circle(float theta) {
    return OpenMesh::Vec2f(cosf(theta), sinf(theta));
}

// Get the initial parameterization of the mesh via
// minimizing the energy E = 1/2 \sum{|| p_i - p_j ||} for the inner point while the boundary points are projected onto the unit circle (Floater Method)
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

    // compute the boundary points
    if (boundaries.size() > 0) {
        std::vector<bool> visited;
        visited.resize(vertices_sz);
        std::fill(visited.begin(), visited.end(), 0);
        double angle = 0.0f;
        double length = 0.0f;
        // compute the total length
        for (auto st = *boundaries.begin(); !visited[st.idx()];) {
            visited[st.idx()] = true;
            decltype(mesh.vv_begin(st)) nx;
            for (nx = mesh.vv_begin(st); nx.is_valid(); ++nx){ 
                if (boundaries.count(*nx) && !visited[(*nx).idx()]) {
                    length += (mesh.point(*nx) - mesh.point(st)).norm();
                    st = *nx;
                    break;
                }
            }
            if (visited[st.idx()]) {
                auto start = boundaries.begin();
                auto i1 = mesh.point(st), i2 = mesh.point(*start);
                length += (i1 - i2).norm();
            }
        }

        // compute the boundary point
        std::fill(visited.begin(), visited.end(), 0);
        res[(*boundaries.begin()).idx()] = project2Circle(0.0f);
        for (auto st = *boundaries.begin(); !visited[st.idx()];) {
            visited[st.idx()] = true;
            for (auto nx = mesh.vv_begin(st); nx.is_valid(); ++nx){ 
                if (boundaries.count(*nx) && !visited[(*nx).idx()]) {
                    angle += 2.0 * M_PI * (mesh.point(*nx) - mesh.point(st)).norm() / length;
                    st = *nx;
                    res[(*nx).idx()] = project2Circle(angle);
                    break;
                }
            }
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
                    sparse.coeffRef(coord * ii,  coord * ii) += 1.0 / dis;
                    sparse.coeffRef(coord * ii + 1,  coord * ii + 1) += 1.0 / dis;
                    B[coord * ii] += p[0] / dis;
                    B[coord * ii + 1] += p[1] / dis;
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
                            sparse.coeffRef(coord * ii + cnt, coord * ij + cnt) -= 0.5 / dis;
                            sparse.coeffRef(coord * ij + cnt, coord * ii + cnt) -= 0.5 / dis;
                            sparse.coeffRef(coord * ii + cnt, coord * ii + cnt) += 1.0 / dis;
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

// Get the energy E_D + E_B for current parameterization and ε
float getEnergy(const Mesh& mesh, const Eigen::VectorXf& params, float epsilon) {
    float energy = 0.0f;
    using OpenMesh::dot;
    // distortion term
    for (auto f: mesh.all_faces()) {
        int i = 0;
        Mesh::VertexHandle vs[3];
        for (auto v: mesh.fv_range(f))
            vs[i++] = v;
        OpenMesh::Vec3f u0(params[2 * vs[0].idx()], params[2 * vs[0].idx() + 1], 0.0f);
        OpenMesh::Vec3f ue1 = OpenMesh::Vec3f(params[2 * vs[1].idx()], params[2 * vs[1].idx() + 1], 0.0f) - u0;
        OpenMesh::Vec3f ue2 = OpenMesh::Vec3f(params[2 * vs[2].idx()], params[2 * vs[2].idx() + 1], 0.0f) - u0;

        OpenMesh::Vec3f p0 = mesh.point(vs[0]);
        OpenMesh::Vec3f pe1 = mesh.point(vs[1]) - p0;
        OpenMesh::Vec3f pe2 = mesh.point(vs[2]) - p0;

        // calculate the area of corresponding triangle
        float su = 0.5 * fabs(OpenMesh::cross(ue1, ue2)[2]);
        float sp = 0.5 * fabs(OpenMesh::cross(pe1, pe2).norm());

        energy += (1 + (sp * sp) / (su * su)) * (sq(ue2.norm()) * sq(pe1.norm()) + sq(ue1.norm()) * sq(pe2.norm()) - 2.0f * dot(ue2, ue1) * dot(pe2, pe1)) / (4.0f * sp);
    }

    // the barrier term
    float tseng = 0.0f;
    for (auto e: mesh.all_edges())
        if (mesh.is_boundary(e)) {
            float seng = 0.0f;
            int st = mesh.to_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
            int end = mesh.from_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
            OpenMesh::Vec2f pst(params[2 * st], params[2 * st + 1]); 
            OpenMesh::Vec2f pend(params[2 * end], params[2 * end + 1]); 
            OpenMesh::Vec2f edge = pend - pst;

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
                            OpenMesh::Vec2f point(params[itr->second * 2], params[itr->second * 2 + 1]);
                            float projected = OpenMesh::dot(point - pst, edge) / edge.norm();
                            float dis = sqrtf(OpenMesh::dot(point - pst, point - pst) - projected * projected);
                            // std::cout << "Dis: " << dis << std::endl;
                            if (dis == 0.0f)
                                return std::numeric_limits<float>::infinity();
                            seng = std::max(seng, epsilon / dis - 1.0f);
                        }
                }
            tseng += seng * seng;
        }
    std::cout << "Rport: Boundary Barrier: [" << tseng << "]" << std::endl <<
        "Distortion: [" << energy << "]\n";
    energy += tseng;
    return energy;
}


// calculate the gradient of the energy function
auto calcGrad(const Mesh& mesh, const Eigen::VectorXf& x, float epsilon) {
    // the energy term
    Eigen::VectorXf gradx(x.size());
    gradx.setZero();
    for (int i = 0; i < gradx.size(); ++i)
        gradx[i] = 0.0f;
    // the gradient for distortion term
    for (auto f: mesh.all_faces()) {
        Mesh::VertexHandle vs[3];
        int i = 0;
        for (auto v: mesh.fv_range(f))
            vs[i++] = v;
        OpenMesh::Vec3f u0(x[2 * vs[0].idx()], x[2 * vs[0].idx() + 1], 0.0f);
        OpenMesh::Vec3f u1(x[2 * vs[1].idx()], x[2 * vs[1].idx() + 1], 0.0f);
        OpenMesh::Vec3f u2(x[2 * vs[2].idx()], x[2 * vs[2].idx() + 1], 0.0f);
        OpenMesh::Vec3f ue1 = u1 - u0;
        OpenMesh::Vec3f ue2 = u2 - u0;

        OpenMesh::Vec3f p0 = mesh.point(vs[0]);
        OpenMesh::Vec3f p1 = mesh.point(vs[1]);
        OpenMesh::Vec3f p2 = mesh.point(vs[2]);
        OpenMesh::Vec3f pe1 = p1 - p0;
        OpenMesh::Vec3f pe2 = p2 - p0;

        // calculate the area of corresponding triangle
        float su = 0.5 * fabs(OpenMesh::cross(ue1, ue2)[2]);
        float sp = 0.5 * fabs(OpenMesh::cross(pe1, pe2).norm());

        float cot0 = OpenMesh::dot(p1 - p0, p2 - p0) / OpenMesh::cross(p1 - p0, p2 - p0).norm();
        float cot1 = OpenMesh::dot(p2 - p1, p0 - p1) / OpenMesh::cross(p2 - p1, p0 - p1).norm();
        float cot2 = OpenMesh::dot(p1 - p2, p0 - p2) / OpenMesh::cross(p1 - p2, p0 - p2).norm();
        
        float term1 = (1 + (sp * sp) / (su * su));
        float term2 = (sq(ue2.norm()) * sq(pe1.norm()) + sq(ue1.norm()) * sq(pe2.norm()) - 2.0f * dot(ue2, ue1) * dot(pe2, pe1)) / (4.0f * sp);

        // dE / du0
        OpenMesh::Vec2f dterm1_du0 = - (sp* sp) / (su * su * su) * OpenMesh::Vec2f(-(u2 - u1)[1], (u2 - u1)[0]);
        OpenMesh::Vec3f dterm2_du0 = -cot2 * u1 - cot1 * u2 + (cot1 + cot2) * u0;
        OpenMesh::Vec2f dE_du0 = dterm1_du0 * term2 + OpenMesh::Vec2f(dterm2_du0[0], dterm2_du0[1]) * term1; 
        gradx[vs[0].idx() * 2] += dE_du0[0];
        gradx[vs[0].idx() * 2 + 1] += dE_du0[1];

        // dE / du1
        OpenMesh::Vec2f dterm1_du1 = - (sp* sp) / (su * su * su) * OpenMesh::Vec2f(-(u0 - u2)[1], (u0 - u2)[0]);
        OpenMesh::Vec3f dterm2_du1 = -cot2 * u0 - cot0 * u2 + (cot0 + cot2) * u1;
        OpenMesh::Vec2f dE_du1 = dterm1_du0 * term2 + OpenMesh::Vec2f(dterm2_du0[0], dterm2_du0[1]) * term1; 
        gradx[vs[1].idx() * 2] += dE_du1[0];
        gradx[vs[1].idx() * 2 + 1] += dE_du1[1];

        // dE / du2
        OpenMesh::Vec2f dterm1_du2 = - (sp* sp) / (su * su * su) * OpenMesh::Vec2f(-(u1 - u0)[1], (u1 - u0)[0]);
        OpenMesh::Vec3f dterm2_du2 = -cot0 * u1 - cot1 * u0 + (cot0 + cot1) * u2;
        OpenMesh::Vec2f dE_du2 = dterm1_du0 * term2 + OpenMesh::Vec2f(dterm2_du0[0], dterm2_du0[1]) * term1; 
        gradx[vs[2].idx() * 2] += dE_du2[0];
        gradx[vs[2].idx() * 2 + 1] += dE_du2[1];
    }
    // the singularity term
    for (auto e: mesh.all_edges())
        if (mesh.is_boundary(e)) {
            float seng = 0.0f;
            int ui_idx = -1;
            float mdis = 0.0f;

            int st = mesh.to_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
            int end = mesh.from_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
            OpenMesh::Vec2f pst(x[2 * st], x[2 * st + 1]);
            OpenMesh::Vec2f pend(x[2 * end], x[2 * end + 1]);
            OpenMesh::Vec2f edge = pend - pst;

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
                            OpenMesh::Vec2f point(x[itr->second * 2], x[itr->second * 2 + 1]);
                            float projected = OpenMesh::dot(point - pst, edge) / edge.norm();
                            float dis = sqrtf(OpenMesh::dot(point - pst, point - pst) - projected * projected);
                            if (dis == 0.0f) {
                                gradx[0] = std::numeric_limits<float>::infinity();
                                return std::move(gradx);
                            }
                            float temp = epsilon / dis - 1.0f;
                            if (temp > seng) {
                                ui_idx = itr->second; 
                                mdis = dis;
                            }
                            seng = std::max(seng, temp);
                        }
                }
            // calculate the gradient
            if (ui_idx != -1 && seng != 0.0f) {
                OpenMesh::Vec2f grad_maxui;
                OpenMesh::Vec2f grad_maxu1;
                OpenMesh::Vec2f grad_maxu2;
                auto u1 = pst;
                auto u2 = pend;
                auto ui = OpenMesh::Vec2f(x[ui_idx * 2], x[ui_idx * 2 + 1]);
                float projected = OpenMesh::dot(ui - u1, u2 - u1) / edge.norm();
                float dis = sqrtf(OpenMesh::dot(ui - u1, ui - u1) - projected * projected);
                // the unit vector othorgonal to the edge
                auto normal = ui - u1 - (u2 - u1).normalize() * projected;
                normal = normal.normalize();
                // dE / du2 u2 is the end point of this edge
                grad_maxu2 = -2 * seng * epsilon * projected / OpenMesh::dot(u2 - u1, u2 - u1) * normal / (mdis * mdis);
                gradx[end * 2] += grad_maxu2[0];
                gradx[end * 2 + 1] += grad_maxu2[1];
                // dE / du1 u1 is the start point of this edge
                projected = OpenMesh::dot(ui - u2, u1 - u2) / edge.norm();
                grad_maxu1 = -2 * seng * epsilon * projected / OpenMesh::dot(u2 - u1, u2 - u1) * normal / (mdis * mdis);
                gradx[st * 2] += grad_maxu1[0];
                gradx[st * 2 + 1] += grad_maxu1[1];
                // dE / dui ui is the point with maximal distance
                grad_maxui = -2 * seng * epsilon * normal / (mdis * mdis);
                gradx[ui_idx * 2] += grad_maxui[0];
                gradx[ui_idx * 2 + 1] += grad_maxui[1];
            }
        }
    return std::move(gradx);
}

// compute the average length of edges
float getAveLen(const Mesh& mesh, const Eigen::VectorXf& x) {
    int count = 0;
    float len = 0.0f;
    for (auto e: mesh.all_edges()) {
        ++count;
        int st = mesh.to_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
        int end = mesh.from_vertex_handle(mesh.halfedge_handle(e, 0)).idx();
        OpenMesh::Vec2f ev(x[2 * st] - x[2 *end], x[2 * st + 1] - x[2 * end + 1]);
        len += ev.norm();
    }

    return len / count; 
}

// initialize the hash map for neighbor search
void initHashMap(const Mesh& mesh, const Eigen::VectorXf& x, float epsilon) {
    neighbors.clear();
    for (int i = 0; i < x.size(); i += 2) {
        OpenMesh::VertexHandle handle(i / 2);
        if (mesh.is_boundary(handle)) {
            OpenMesh::Vec2i idx(x[i] / epsilon, x[i + 1] / epsilon);
            neighbors.insert(std::make_pair(idx, i / 2));
        }
    }
}

// get the step length for line search
double getStepLength(const Mesh& mesh, const Eigen::VectorXf& params, const Eigen::VectorXf& p, const Eigen::VectorXf& grad, const float& epsilon) {
    double alpha0 = 1.0f;
    constexpr double c1 = 1e-4;
    double ene0 = getEnergy(mesh, params, epsilon);
    double enei = getEnergy(mesh, params + alpha0 * p, epsilon);
    double dphi = p.dot(grad);
    while (enei > ene0) {
        alpha0 /= 2.0;
        enei = getEnergy(mesh, params + alpha0 * p, epsilon);
    }
    double inite = enei, inite0 = ene0;
    double init = alpha0;
    if (enei < ene0 + c1 * alpha0 * dphi)
        return  alpha0;
    double alpha = - dphi * alpha0 * alpha0 / (2.0f * (enei - ene0 - dphi * alpha0));
    double enei_0 = enei; // energy_{i-1}
    enei = getEnergy(mesh, params + alpha * p, epsilon);
    std::cout << alpha << " " << enei << std::endl;
    while (!(enei < ene0 + c1 * alpha * dphi)) {
        // calculate the coeff of the cubic interpolation
        double dim0 = enei - ene0 - dphi * alpha;
        double dim1 = enei_0 - ene0 - dphi * alpha0;
        double numerator = 1.0f / (alpha0 * alpha0 * alpha * alpha * (alpha - alpha0));

        double a = (alpha0 * alpha0 * dim0 - alpha * alpha * dim1) * numerator;
        double b = (-alpha0 * alpha0 * alpha0 * dim0 + alpha * alpha * alpha * dim1) * numerator;
        // calculate the aplha
        alpha0 = alpha;
        alpha = -b + sqrtf(b * b - 3 * a *dphi) / (3 * a);
        // update 
        enei_0 = enei;
        enei = getEnergy(mesh, params + alpha * p, epsilon);
    }
    return alpha;
}

// To be tested
auto LBFGS(Mesh& mesh, Eigen::VectorXf x, const double& error) {
    static constexpr int sz = 7;
    auto test = [&error](const Eigen::VectorXf& x) -> bool {
        return x.norm() > error;
    };

    Eigen::VectorXf ss[sz], ys[sz];
    float rhos[sz], alphas[sz];
    float epsilon = getAveLen(mesh, x) * 0.25f;
    Eigen::VectorXf gradx = calcGrad(mesh, x, epsilon);
    int k = 0;
    do{
        // compute the average length
        epsilon = getAveLen(mesh, x) * 0.25f;
        // initialize the hash map for neighborhood query
        initHashMap(mesh, x, epsilon);
        auto q = gradx;
        auto gradx0 = gradx;
        Eigen::VectorXf r(x.size());
        for (int i = k - 1; i >= std::max(0, k - sz); --i) {
            alphas[i % sz] = q.dot(ss[i % sz]) * rhos[i % sz];
            q -= alphas[i % sz] * ys[i % sz]; 
        }
        if (k == 0)
            r = q;
        else r = ss[(k - 1) % sz].dot(ys[(k - 1) % sz]) / ys[(k - 1) % sz].dot(ys[(k - 1) % sz]) * q;
        for (int i = std::max(0, k - sz); i < k; ++i) {
            auto beta = rhos[i % sz] * r.dot(ys[i %sz]);
            r += ss[i] * (alphas[i % sz] - beta);
        }
        float step_len = getStepLength(mesh, x, -r, gradx0, epsilon);
        ss[k % sz] = step_len * -r;
        x += ss[k % sz];
        gradx = calcGrad(mesh, x, epsilon);
        ys[k % sz] = gradx - gradx0;
        rhos[k % sz] = 1.0f / ss[k % sz].dot(ys[k % sz]);
        ++k;
    }while(test(gradx)); // todo!

    return std::move(x);
}

std::vector<OpenMesh::Vec2f> paratimization(Mesh& mesh, double error) {
    auto init_param = init(mesh); // the initial parameterization for the mesh
    Eigen::VectorXf x(init_param.size() * 2);
    for (int i = 0; i < init_param.size(); ++i) {
        x[2 * i] = init_param[i][0];
        x[2 * i + 1] = init_param[i][1];
    }

    LBFGS(mesh, x, error); // optimize the parameterization using Quasi Newton Method

    for (int i = 0; i < init_param.size(); ++i) {
        init_param[i][0] = x[2 * i];
        init_param[i][1] = x[2 * i + 1];
    }

    return std::move(init_param);
}
