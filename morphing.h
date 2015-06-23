#ifndef MORPHING_H
#define MORPHING_H

#include <Eigen/Core>

struct MorphParams
{
    MorphParams();

    int framesPerMorph;
    int DFTCoeffs;
    double samplingResolution;
    double reconstructionResolution;
    double regularization;
};

class Morphing
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Morphing(MorphParams params, const Eigen::VectorXd &polygon1, const Eigen::VectorXd &polygon2);

    void renderInterpolation(int tick);

private:
    double area(const Eigen::VectorXd &pts);
    Eigen::Vector2d centroid(const Eigen::VectorXd &pts);
    Eigen::Matrix2d principalAxes(const Eigen::VectorXd &pts);
    void subsample(const Eigen::VectorXd &inpts, const Eigen::Vector2d &centroid, Eigen::VectorXd &outpts);
    void toPolar(const Eigen::VectorXd &inpts, Eigen::VectorXd &outpts);
    void slowDFT(const Eigen::VectorXd &inputpts, Eigen::VectorXd &coeffs);
    void inverseDFTMatrix(Eigen::MatrixXd &M);
    void inverseDFT(const Eigen::VectorXd &coeffs, const Eigen::Vector2d &centroid, Eigen::VectorXd &pts);

    MorphParams params_;
    Eigen::Vector2d centroid_;
    Eigen::Vector2d translation_;

    Eigen::VectorXd coeffs1_;
    Eigen::VectorXd coeffs2_;

    Eigen::MatrixXd M_;
};

#endif // MORPHING_H
