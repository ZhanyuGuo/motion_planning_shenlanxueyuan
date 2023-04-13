#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;
using namespace Eigen;

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint() {}
TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint() {}

// define factorial function, input i, output i!
int TrajectoryGeneratorWaypoint::Factorial(int x)
{
    int fac = 1;
    for (int i = x; i > 0; i--)
        fac = fac * i;

    return fac;
}

void TrajectoryGeneratorWaypoint::GetHessian(const int n_seg,
                                             const int d_order,
                                             const Eigen::VectorXd &Time,
                                             Eigen::SparseMatrix<double> &hession)
{
    int p_order = 2 * d_order - 1;
    int p_num1d = p_order + 1;
    hession.resize(n_seg * p_num1d, n_seg * p_num1d);
    hession.setZero();

    for (int k = 0; k < n_seg; ++k)
    {
        for (int i = d_order; i < p_num1d; ++i)
        {
            for (int j = d_order; j < p_num1d; ++j)
            {
                double value = 1.0 * Factorial(i) / Factorial(i - d_order) * Factorial(j) / Factorial(j - d_order) / (i + j - 2 * d_order + 1) * pow(Time(k), i + j - 2 * d_order + 1);
                hession.insert(k * p_num1d + i, k * p_num1d + j) = value;
            }
        }
    }
}

void TrajectoryGeneratorWaypoint::InsertCoff(const int row,
                                             const int col,
                                             Eigen::SparseMatrix<double> &linearMatrix,
                                             const double t,
                                             const int d_order,
                                             bool one_line,
                                             bool reverse)
{
    int p_num1d = 2 * d_order;

    int flag = d_order;
    if (one_line)
        flag = 1;

    Eigen::MatrixXd coff(d_order, p_num1d);

    if (d_order == 4)
    {
        coff << 1.0, 1.0 * t, 1.0 * pow(t, 2), 1.0 * pow(t, 3), 1.0 * pow(t, 4), 1.0 * pow(t, 5), 1.0 * pow(t, 6), 1.0 * pow(t, 7),
            0.0, 1.0, 2.0 * t, 3.0 * pow(t, 2), 4.0 * pow(t, 3), 5.0 * pow(t, 4), 6.0 * pow(t, 5), 7.0 * pow(t, 6),
            0.0, 0.0, 2.0, 6.0 * t, 12.0 * pow(t, 2), 20.0 * pow(t, 3), 30.0 * pow(t, 4), 42.0 * pow(t, 5),
            0.0, 0.0, 0.0, 6.0, 24.0 * t, 60.0 * pow(t, 2), 120.0 * pow(t, 3), 210.0 * pow(t, 4);
    }
    else if (d_order == 3)
    {
        coff << 1.0, 1.0 * t, 1.0 * pow(t, 2), 1.0 * pow(t, 3), 1.0 * pow(t, 4), 1.0 * pow(t, 5),
            0.0, 1.0, 2.0 * t, 3.0 * pow(t, 2), 4.0 * pow(t, 3), 5.0 * pow(t, 4),
            0.0, 0.0, 2.0, 6.0 * t, 12.0 * pow(t, 2), 20.0 * pow(t, 3);
    }
    else
    {
        ROS_ERROR("Only d_order = 3 or 4 supported!");
    }

    if (reverse)
        coff = coff * (-1.0);

    for (int i = 0; i < d_order && i < flag; ++i)
    {
        for (int j = 0; j < p_num1d; ++j)
        {
            linearMatrix.insert(row + i, col + j) = coff(i, j);
        }
    }
}

void TrajectoryGeneratorWaypoint::GetLinearConstraintsMatrix(const int n_seg,
                                                             const int d_order,
                                                             const Eigen::VectorXd &Time,
                                                             Eigen::SparseMatrix<double> &linearMatrix)
{
    int p_order = 2 * d_order - 1;
    int p_num1d = p_order + 1;
    linearMatrix.resize(2 * d_order + (n_seg - 1) * (d_order + 1), p_num1d * n_seg);

    // start and end condition
    int row = 0, col = 0;
    InsertCoff(row, col, linearMatrix, 0, d_order, false, false);

    row += d_order;
    col = (n_seg - 1) * p_num1d;
    InsertCoff(row, col, linearMatrix, Time(n_seg - 1), d_order, false, false);

    // waypoint positions
    row += d_order;
    for (int k = 0; k < n_seg - 1; k++)
    {
        InsertCoff(row + k, k * p_num1d, linearMatrix, Time(k), d_order, true, false);
    }

    // continuty
    row += n_seg - 1;
    for (int k = 0; k < n_seg - 1; k++)
    {
        InsertCoff(row, k * p_num1d, linearMatrix, Time(k), d_order, false, false);
        InsertCoff(row, (k + 1) * p_num1d, linearMatrix, 0, d_order, false, true);
        row += d_order;
    }
}

/**
 *  STEP 2: Learn the "Closed-form solution to minimum snap" in L5, then finish this PolyQPGeneration function
 */
Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
    const int d_order,           // the order of derivative
    const Eigen::MatrixXd &Path, // waypoints coordinates (3d)
    const Eigen::MatrixXd &Vel,  // boundary velocity
    const Eigen::MatrixXd &Acc,  // boundary acceleration
    const Eigen::VectorXd &Time) // time allocation in each segment
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order = 2 * d_order - 1; // the order of polynomial
    int p_num1d = p_order + 1;     // the number of variables in each segment

    int m = Time.size();                                 // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d); // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);

    // calculate Q
    MatrixXd Q = MatrixXd::Zero(p_num1d * m, p_num1d * m);
    for (int k = 0; k < m; k++)
    {
        MatrixXd Q_k = MatrixXd::Zero(p_num1d, p_num1d);
        for (int i = d_order; i < p_num1d; i++)
        {
            for (int j = d_order; j < p_num1d; j++)
            {
                // Q_j(i + 1, k + 1) = factorial(i) / factorial(i - 4) * factorial(k) / factorial(k - 4) / (i + k - n_order) * ts(j) ^ (i + k - n_order);
                Q_k(i, j) = 1.0 * Factorial(i) / Factorial(i - d_order) * Factorial(j) / Factorial(j - d_order) / (i + j - p_order) * pow(Time(k), i + j - p_order);
            }
        }
        Q.block(k * p_num1d, k * p_num1d, p_num1d, p_num1d) = Q_k;
    }

    // calculate M
    MatrixXd M = MatrixXd::Zero(p_num1d * m, p_num1d * m);
    MatrixXd coeff(d_order, p_num1d);
    coeff << 1, 1, 1, 1, 1, 1, 1, 1,
        0, 1, 2, 3, 4, 5, 6, 7,
        0, 0, 2, 6, 12, 20, 30, 42,
        0, 0, 0, 6, 24, 60, 120, 210;

    for (int k = 0; k < m; k++)
    {
        MatrixXd M_k = MatrixXd::Zero(p_num1d, p_num1d);
        double t = Time(k);
        for (int i = 0; i < d_order; i++)
        {
            M_k(i, i) = coeff(i, i);
        }

        for (int i = 0; i < d_order; i++)
        {
            for (int j = i; j < p_num1d; j++)
            {
                if (i == j)
                {
                    M_k(i + d_order, j) = coeff(i, j);
                }
                else
                {
                    M_k(i + d_order, j) = coeff(i, j) * pow(t, j - i);
                }
            }
        }
        M.block(k * p_num1d, k * p_num1d, p_num1d, p_num1d) = M_k;
    }

    // calculate Ct
    int ct_rows = 2 * d_order * m;
    int ct_cols = 2 * d_order + (m - 1) * d_order;
    MatrixXd Ct = MatrixXd::Zero(ct_rows, ct_cols);
    // hash map
    vector<int> d_vector;
    for (int k = 0; k < m; k++)
    {
        for (int t = 0; t < 2; t++)
        {
            for (int d = 0; d < d_order; d++)
            {
                d_vector.push_back(k * 100 + t * 10 + d);
            }
        }
    }

    int val, row;
    int k, t, d;
    int col = 0;

    // start condition
    k = 0;
    t = 0;
    for (d = 0; d < d_order; d++)
    {
        val = k * 100 + t * 10 + d;
        auto it = find(d_vector.begin(), d_vector.end(), val);
        row = distance(d_vector.begin(), it);
        Ct(row, col) = 1;
        col++;
    }

    // waypoint position
    t = 1;
    d = 0;
    for (k = 0; k < m - 1; k++)
    {
        val = k * 100 + t * 10 + d;
        auto it = find(d_vector.begin(), d_vector.end(), val);
        row = distance(d_vector.begin(), it);
        Ct(row, col) = 1;

        val = (k + 1) * 100 + (t - 1) * 10 + d;
        it = find(d_vector.begin(), d_vector.end(), val);
        row = distance(d_vector.begin(), it);
        Ct(row, col) = 1;

        col++;
    }

    // end condition
    k = m - 1;
    t = 1;
    for (d = 0; d < d_order; d++)
    {
        val = k * 100 + t * 10 + d;
        auto it = find(d_vector.begin(), d_vector.end(), val);
        row = distance(d_vector.begin(), it);
        Ct(row, col) = 1;
        col++;
    }

    // other free, but need continuty
    t = 1;
    for (k = 0; k < m - 1; k++)
    {
        for (d = 1; d < d_order; d++)
        {
            val = k * 100 + t * 10 + d;
            auto it = find(d_vector.begin(), d_vector.end(), val);
            row = distance(d_vector.begin(), it);
            Ct(row, col) = 1;

            val = (k + 1) * 100 + (t - 1) * 10 + d;
            it = find(d_vector.begin(), d_vector.end(), val);
            row = distance(d_vector.begin(), it);
            Ct(row, col) = 1;

            col++;
        }
    }

    MatrixXd C = Ct.transpose();
    MatrixXd M_inv = M.inverse();
    MatrixXd M_inv_t = M_inv.transpose();
    MatrixXd R = C * M_inv_t * Q * M_inv * Ct;

    int num_d_F = 2 * d_order + m - 1;
    int num_d_P = (m - 1) * (d_order - 1);

    MatrixXd R_FP = R.topRightCorner(num_d_F, num_d_P);
    MatrixXd R_PP = R.bottomRightCorner(num_d_P, num_d_P);

    for (int dim = 0; dim < 3; dim++)
    {
        VectorXd wayPoints = Path.col(dim);
        VectorXd d_F = VectorXd::Zero(num_d_F);

        d_F(0) = wayPoints(0);
        for (int i = 0; i < m - 1; i++)
        {
            d_F(i + d_order) = wayPoints(i + 1);
        }
        d_F(m - 1 + d_order) = wayPoints(m);

        VectorXd d_P = -1.0 * R_PP.inverse() * R_FP.transpose() * d_F;
        VectorXd d_total(num_d_F + num_d_P);
        d_total << d_F, d_P;

        VectorXd poly_coef_1d = M.inverse() * Ct * d_total;
        MatrixXd poly_coef_1d_t = poly_coef_1d.transpose();

        for (int k = 0; k < m; k++)
        {
            PolyCoeff.block(k, dim * p_num1d, 1, p_num1d) = poly_coef_1d_t.block(0, k * p_num1d, 1, p_num1d);
        }
    }

    return PolyCoeff;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
    const int d_order,           // the order of derivative
    const Eigen::MatrixXd &Path, // waypoints coordinates (3d)
    const Eigen::MatrixXd &Vel,  // boundary velocity
    const Eigen::MatrixXd &Acc,  // boundary acceleration
    const Eigen::VectorXd &Time, // time allocation in each segment
    OsqpEigen::Solver &solver)
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order = 2 * d_order - 1; // the order of polynomial
    int p_num1d = p_order + 1;     // the number of variables in each segment

    int m = Time.size();                                 // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d); // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);

    // set number of values and constraints
    solver.data()->setNumberOfVariables(m * p_num1d);
    solver.data()->setNumberOfConstraints(2 * d_order + (m - 1) * (d_order + 1));

    // set Hession
    Eigen::SparseMatrix<double> hessian;
    GetHessian(m, d_order, Time, hessian);
    if (!solver.data()->setHessianMatrix(hessian))
    {
        ROS_ERROR("failed to set Hessian");
        return PolyCoeff;
    }

    // set linear constraints
    Eigen::SparseMatrix<double> linearMatrix;
    GetLinearConstraintsMatrix(m, d_order, Time, linearMatrix);
    if (!solver.data()->setLinearConstraintsMatrix(linearMatrix))
    {
        ROS_ERROR("failed to set linearMatrix");
        return PolyCoeff;
    }

    // set gradient
    Eigen::VectorXd gradient(p_num1d * m);
    gradient.setZero();
    if (!solver.data()->setGradient(gradient))
    {
        ROS_ERROR("failed to set gradient");
        return PolyCoeff;
    }

    // boundary
    Eigen::VectorXd lowbound = VectorXd::Zero(2 * d_order + (m - 1) * (d_order + 1));
    Eigen::VectorXd upbound = VectorXd::Zero(2 * d_order + (m - 1) * (d_order + 1));

    solver.data()->setLowerBound(lowbound);
    solver.data()->setUpperBound(upbound);

    if (!solver.isInitialized())
    {
        solver.initSolver();
    }

    for (int dim = 0; dim < 3; dim++)
    {
        VectorXd wayPoints = Path.col(dim);

        // start condition
        lowbound(0) = wayPoints(0);
        upbound(0) = wayPoints(0);

        // end condition
        lowbound(d_order) = wayPoints(m);
        upbound(d_order) = wayPoints(m);

        // waypoint position
        for (int i = 0; i < m - 1; i++)
        {
            lowbound(2 * d_order + i) = wayPoints(i + 1);
            upbound(2 * d_order + i) = wayPoints(i + 1);
        }

        solver.updateBounds(lowbound, upbound);
        solver.solveProblem();
        Eigen::VectorXd poly_coef_1d = solver.getSolution();
        MatrixXd poly_coef_1d_t = poly_coef_1d.transpose();

        for (int k = 0; k < m; k++)
        {
            PolyCoeff.block(k, dim * p_num1d, 1, p_num1d) = poly_coef_1d_t.block(0, k * p_num1d, 1, p_num1d);
        }
    }

    // reset
    solver.data()->clearHessianMatrix();
    solver.data()->clearLinearConstraintsMatrix();
    solver.clearSolverVariables();
    solver.clearSolver();

    return PolyCoeff;
}