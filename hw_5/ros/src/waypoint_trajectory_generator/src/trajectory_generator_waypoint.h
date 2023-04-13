#ifndef _TRAJECTORY_GENERATOR_WAYPOINT_H_
#define _TRAJECTORY_GENERATOR_WAYPOINT_H_

#include <Eigen/Eigen>
#include <vector>
#include <OsqpEigen/OsqpEigen.h>

class TrajectoryGeneratorWaypoint
{
private:
    double _qp_cost;
    Eigen::MatrixXd _Q;
    Eigen::VectorXd _Px, _Py, _Pz;

public:
    TrajectoryGeneratorWaypoint();
    ~TrajectoryGeneratorWaypoint();

    Eigen::MatrixXd PolyQPGeneration(
        const int order,
        const Eigen::MatrixXd &Path,
        const Eigen::MatrixXd &Vel,
        const Eigen::MatrixXd &Acc,
        const Eigen::VectorXd &Time);

    Eigen::MatrixXd PolyQPGeneration(
        const int d_order,
        const Eigen::MatrixXd &Path,
        const Eigen::MatrixXd &Vel,
        const Eigen::MatrixXd &Acc,
        const Eigen::VectorXd &Time,
        OsqpEigen::Solver &solver);

    int Factorial(int x);

    void GetHessian(const int n_seg,
                    const int d_order,
                    const Eigen::VectorXd &Time,
                    Eigen::SparseMatrix<double> &hession);

    void InsertCoff(const int row,
                    const int col,
                    Eigen::SparseMatrix<double> &linearMatrix,
                    const double t,
                    const int d_order,
                    bool one_line,
                    bool reverse);

    void GetLinearConstraintsMatrix(const int n_seg,
                                    const int d_order,
                                    const Eigen::VectorXd &Time,
                                    Eigen::SparseMatrix<double> &linearMatrix);
};

#endif
