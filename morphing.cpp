#include "morphing.h"
#include <QGLWidget>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

const double PI = 3.1415926535898;

MorphParams::MorphParams()
{
    samplingResolution = 5000;
    framesPerMorph = 100;
    DFTCoeffs = 1000;
    reconstructionResolution = 1000;
    regularization = 1e-3;
}

Morphing::Morphing(MorphParams params, const Eigen::VectorXd &polygon1, const Eigen::VectorXd &polygon2) : params_(params)
{
    Vector2d centroid1 = centroid(polygon1);
    Vector2d centroid2 = centroid(polygon2);
    centroid_ = centroid1;
    translation_ = centroid2 - centroid1;

    VectorXd samples1, samples2;
    subsample(polygon1, centroid1, samples1);
    subsample(polygon2, centroid2, samples2);

    VectorXd polar1, polar2;
    toPolar(samples1, polar1);
    toPolar(samples2, polar2);

    slowDFT(polar1, coeffs1_);
    slowDFT(polar2, coeffs2_);

    inverseDFTMatrix(M_);
}

double Morphing::area(const VectorXd &pts)
{
    int npts = (int)pts.size()/2;

    double tot = 0;

    for(int i=0; i<npts; i++)
    {
        int next = (i+1)%npts;
        Vector2d v1 = pts.segment<2>(2*i);
        Vector2d v2 = pts.segment<2>(2*next) - v1;

        double term = 0.5*(-v1[0]*v2[1]+v1[1]*v2[0]);
        tot += term;
    }
    return tot;
}

Vector2d Morphing::centroid(const VectorXd &pts)
{
    Vector2d result(0,0);
    int npts = (int)pts.size()/2;
    for(int i=0; i<npts; i++)
    {
        int next = (i+1)%npts;
        Vector2d p1 = pts.segment<2>(2*i);
        Vector2d p2 = pts.segment<2>(2*next);
        result[0] += 0.5 * (p1[0]*p1[0]/3.0 + p1[0]*p2[0]/3.0 + p2[0]*p2[0]/3.0)*(p1[1] - p2[1]);
        result[1] += 0.5 * (p1[1]*p1[1]/3.0 + p1[1]*p2[1]/3.0 + p2[1]*p2[1]/3.0)*(p2[0] - p1[0]);
    }
    result /= area(pts);
    return result;
}

void Morphing::renderInterpolation(int tick)
{
    double t = double(tick % params_.framesPerMorph)/double(params_.framesPerMorph - 1);
    VectorXd polygon;
    VectorXd interpcoeffs = ((1-t)*coeffs1_+t*coeffs2_);///sqrt(t*t+(1-t)*(1-t));
    inverseDFT(interpcoeffs, centroid_, polygon);

    int npts = (int)polygon.size()/2;
    glColor3f(0.0, 0.0, 1.0);
    glBegin(GL_LINE_STRIP);

    for(int i=0; i<npts; i++)
    {
        glVertex2d(polygon[2*i] + t*translation_[0], polygon[2*i+1] + t*translation_[1]);
    }
    if(npts > 0)
        glVertex2d(polygon[0] + t*translation_[0], polygon[1] + t*translation_[1]);
    glEnd();
}

void Morphing::subsample(const VectorXd &inpts, const Vector2d &centroid, VectorXd &outpts)
{
    std::vector<int> samplesperseg;
    int totsamples = 0;

    double arclength = 0;

    int numpts = (int)inpts.size()/2;
    for(int i=0; i<numpts; i++)
    {
        int next = (i+1)%numpts;
        Vector2d seg = inpts.segment<2>(2*next)-inpts.segment<2>(2*i);
        double len = seg.norm();
        arclength += len;
    }

    for(int i=0; i<numpts; i++)
    {
        int next = (i+1)%numpts;
        Vector2d seg = inpts.segment<2>(2*next)-inpts.segment<2>(2*i);
        double len = seg.norm();
        int samples = ceil(len/arclength*double(params_.samplingResolution));
        samplesperseg.push_back(samples);
        totsamples += samples;
    }

    outpts.resize(2*totsamples);

    int cur = 0;

    for(int i=0; i<numpts; i++)
    {
        int samples = samplesperseg[i];
        int next = (i+1)%numpts;
        Vector2d pt1 = inpts.segment<2>(2*i);
        Vector2d seg = inpts.segment<2>(2*next)-inpts.segment<2>(2*i);
        double step = 1.0/samples;
        for(int j=0; j<samples; j++)
        {
            outpts[cur++] = pt1[0] + j*step*seg[0] - centroid[0];
            outpts[cur++] = pt1[1] + j*step*seg[1] - centroid[1];
        }
    }    
    assert(cur == 2*totsamples);
}

void Morphing::toPolar(const VectorXd &inpts, VectorXd &outpts)
{
    int npts = (int)inpts.size()/2;
    outpts.resize(2*npts);

    for(int i=0; i<npts; i++)
    {
        Vector2d pt = inpts.segment<2>(2*i);
        outpts[2*i] = pt.norm();
        outpts[2*i+1] = atan2(pt[1], pt[0]);
    }
}

void Morphing::slowDFT(const VectorXd &inputpts, VectorXd &coeffs)
{
    int ncoeffs = params_.DFTCoeffs;
    int npts = (int)inputpts.size()/2;
    MatrixXd M(npts, ncoeffs);
    VectorXd rhs(npts);
    M.setZero();
    for(int i=0; i<npts; i++)
    {
        for(int j=0; j<ncoeffs; j++)
        {
            double coeff;
            if(j%2 == 0)
            {
                coeff = cos(inputpts[2*i+1]*j/2);
            }
            else
                coeff = sin(inputpts[2*i+1]*(j+1)/2);

            M.coeffRef(i, j) = coeff;
        }
        rhs[i] = inputpts[2*i];
    }

    VectorXd MTrhs = M.transpose()*rhs;
    MatrixXd MTM = M.transpose()*M;
    for(int i=0; i<ncoeffs; i++)
    {
        MTM.coeffRef(i,i) += params_.regularization*i*i/4.0;
    }

    coeffs = LDLT<MatrixXd>(MTM).solve(MTrhs);
}

void Morphing::inverseDFTMatrix(MatrixXd &M)
{
    M.resize(params_.reconstructionResolution, params_.DFTCoeffs);
    M.setZero();
    for(int i=0; i<params_.reconstructionResolution; i++)
    {
        for(int j=0; j<params_.DFTCoeffs; j++)
        {
            double coeff;
            if(j%2 == 0)
                coeff = cos(2*PI*i/params_.reconstructionResolution*j/2);
            else
                coeff = sin(2*PI*i/params_.reconstructionResolution*(j+1)/2);
            M.coeffRef(i,j) = coeff;
        }
    }
}

void Morphing::inverseDFT(const VectorXd &coeffs, const Vector2d &centroid, VectorXd &pts)
{
    pts.resize(2*params_.reconstructionResolution);
    VectorXd rs = M_*coeffs;
    for(int i=0; i<params_.reconstructionResolution; i++)
    {
        double theta = 2*PI*i/params_.reconstructionResolution;
        pts[2*i] = centroid[0] + rs[i]*cos(theta);
        pts[2*i+1] = centroid[1] + rs[i]*sin(theta);
    }
}
