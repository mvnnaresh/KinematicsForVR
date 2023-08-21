#pragma once

#include <iomanip>
#include <kdl/kdl.hpp>
#include <kdl/chain.hpp>
#include <kdl/chainfksolver.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/chainiksolvervel_wdls.hpp>
#include <kdl/chainiksolverpos_nr.hpp>
#include <kdl/chainiksolverpos_nr_jl.hpp>
#include <visp3/core/vpHomogeneousMatrix.h>

#define M_PI 3.14159265358979323846
#define M_PI_2 1.57079632679489661923
#define DLS_LAMBDA 0.3
#define RADIANS(A)  A * (M_PI / 180)
#define DEGREES(A)  A * (180 / M_PI)

struct ChainData
{
    KDL::Chain chain;
    std::vector<double> q_min;
    std::vector<double> q_max;
    std::vector<double> q_vmax;
};

enum IKMethod
{
    KDL_PINV_IK = 0, //fastest of all
    JPINV_NSO,
    DLS_NSO
};

#define IK_METHOD KDL_PINV_IK

inline ChainData robotdata()
{
    ChainData kukaiiwa_14_r820;

    kukaiiwa_14_r820.chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ), KDL::Frame(KDL::Rotation::RotZ(0), KDL::Vector(0.0, 0.0, 0.36))));
    kukaiiwa_14_r820.chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotY), KDL::Frame(KDL::Rotation::RotZ(-M_PI_2))));
    kukaiiwa_14_r820.chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ), KDL::Frame(KDL::Rotation::RotZ(-M_PI_2), KDL::Vector(0.0, 0.0, 0.42))));
    kukaiiwa_14_r820.chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotY), KDL::Frame(KDL::Rotation::RotZ(M_PI_2))));
    kukaiiwa_14_r820.chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ), KDL::Frame(KDL::Rotation::RotZ(M_PI_2), KDL::Vector(0.0, 0.0, 0.4))));
    kukaiiwa_14_r820.chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotY), KDL::Frame(KDL::Rotation::RotZ(M_PI_2))));
    kukaiiwa_14_r820.chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ), KDL::Frame(KDL::Rotation::RotZ(-M_PI_2), KDL::Vector(0.0, 0.0, 0.126))));

    kukaiiwa_14_r820.q_min.assign(
    {
        RADIANS(-170.0), RADIANS(-120.0), RADIANS(-170.0), RADIANS(-120.0), RADIANS(-170.0), RADIANS(-120.0), RADIANS(-175.0)
    });
    kukaiiwa_14_r820.q_max.assign(
    {
        RADIANS(170.0), RADIANS(120.0), RADIANS(170.0), RADIANS(120.0), RADIANS(170.0), RADIANS(120.0), RADIANS(175.0)
    });

    kukaiiwa_14_r820.q_vmax.assign(
    {
        RADIANS(98), RADIANS(98), RADIANS(100), RADIANS(130), RADIANS(140), RADIANS(180), RADIANS(180)
    });

    return kukaiiwa_14_r820;
}

class dxKin
{
public:
    dxKin();

    ~dxKin() {}

    int numJoints{};

    vpHomogeneousMatrix getForwardPosition(std::vector <double> j);

    //Input in DEGREES -  output is in DEGREES
    std::vector <double> getJointCommand(std::vector<double> qcurr, vpHomogeneousMatrix H);

protected:
    std::vector<double> low;

    std::vector<double> high;

    std::vector<double> velLim;

private:
    ChainData robot;

    std::shared_ptr<KDL::JntArray> q;

    std::shared_ptr<KDL::JntArray> qdot;

    std::shared_ptr<KDL::JntArray> qdotdot;

    KDL::ChainFkSolverPos_recursive *fksolver;

    KDL::ChainIkSolverVel_pinv *iksolver1v;

    KDL::ChainIkSolverPos_NR_JL *iksolver;

    void setJointPos(std::vector<double> pos);

    vpMatrix getBaseJacobianMatrix(std::vector<double> q);

    std::vector<double> getPositionWithNullSpace(std::vector<double> xdot, std::vector<double> qcurrent,
            std::vector<double> nullspaceQ0_dot);

    std::vector<double> getPositionWithNullSpaceAndDLS(std::vector<double> xdot, std::vector<double> qcurrent,
            std::vector<double> nullspaceQ0_dot);

    double getManipulability(std::vector<double> qcurrent);

    //note: 'q' values in DEGREES
    KDL::Frame fkSolver(std::vector<double> q);

    //Note: 'seed' configuration and OUTPUT are in DEGREES
    std::vector<double> ikSolverKDL(vpHomogeneousMatrix cartPos, std::vector<double> seed);
    std::vector<double> ikSolverNSO(vpHomogeneousMatrix cartPos, std::vector<double> seed);

    vpColVector computePoseError(vpHomogeneousMatrix Hd, vpHomogeneousMatrix Hc);

    vpColVector computeQuaternionError(vpQuaternionVector Qd, vpQuaternionVector Qc);

    KDL::Frame vpHomMatToKDLFrame(vpHomogeneousMatrix cartPos)
    {
        KDL::Rotation RF = KDL::Rotation(cartPos[0][0], cartPos[0][1], cartPos[0][2],
                                         cartPos[1][0], cartPos[1][1], cartPos[1][2],
                                         cartPos[2][0], cartPos[2][1], cartPos[2][2]);

        KDL::Vector TF = KDL::Vector(cartPos[0][3], cartPos[1][3], cartPos[2][3]);
        KDL::Frame F = KDL::Frame(RF, TF);

        return F;
    }

    vpHomogeneousMatrix KDLFrameTovpHomMat(KDL::Frame F)
    {
        vpHomogeneousMatrix H;

        vpRotationMatrix R;
        vpTranslationVector T;

        for (int i = 0; i < 3; i++)
        {
            for (int k = 0; k < 3; k++)
            {
                R[i][k] = F.M(i, k);
            }
            T[i] = F.p(i);
        }

        H.buildFrom(T, R);

        return H;
    }

    std::vector<double> KDLJntArrayToVector(KDL::JntArray ja)
    {
        std::vector<double> res;
        for (unsigned int i = 0; i < ja.rows(); i++)
        {
            res.push_back(ja(i));
        }

        return res;
    }

    std::vector<double> vectorToRadians(std::vector<double> q)
    {
        std::vector<double>  s;
        s.resize(q.size());

        for (auto i = 0; i < q.size(); i++)
        {
            s[i] = vpMath::rad(q[i]);
        }
        return s;
    }

    std::vector<double> vectorToDegrees(std::vector<double> q)
    {
        std::vector<double>  s;
        s.resize(q.size());

        for (auto i = 0; i < q.size(); i++)
        {
            s[i] = vpMath::deg(q[i]);
        }
        return s;
    }
};

