#pragma once
#include "dxKin.h"


#define IK_METHOD IKMethod::DLS_NSO   //KDL_PINV_IK 0, JPINV_NSO, DLS_NSO

class dxKinInterface
{
public:
    dxKinInterface();

    /************************************************************************
     * INPUT:
     * current joints: q1, ..... q7 in degrees
     * target position: X, Y, Z, Q.x, Q.y, Q.z, Q.w (position + quaternions)
     * position in meters
     *
     * OUTPUT:
     * Vector of joint angles in degrees (if IK solution found).
     * If NO IK solution found, then returns current joint angles
     ************************************************************************/
    std::vector<double> getJointCommand(std::vector<double> currentJoints, const std::vector<double>& targetPosition);

    void test();
protected:
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

private:
    std::shared_ptr<dxKin> kin;

};

