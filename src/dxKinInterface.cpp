#include <utility>

#include "..\include\dxKinInterface.h"

dxKinInterface::dxKinInterface()
{
    kin = std::make_shared<dxKin>(this->robotdata());

    kin->setIKMethod(IK_METHOD);
}

std::vector<double> dxKinInterface::getJointCommand(std::vector<double> currentJoints, const std::vector<double>& targetPosition)
{
    if(static_cast<int>(currentJoints.size()) != kin->numJoints || targetPosition.size() != 7)
    {
        printf("[Error] >> Input vectors size does not satisfy the working criteria.\n\n");
        printf("[Help Message ] >> \n currentJoints -> size should be same as robot degrees of freedom\n targetPosition -> size should be 7");
        return currentJoints;
    }

    else
    {
        vpTranslationVector T(targetPosition[0], targetPosition[1], targetPosition[2]);
        vpQuaternionVector Q(targetPosition[3], targetPosition[4], targetPosition[5], targetPosition[6]);

        vpRotationMatrix R(Q);

        vpHomogeneousMatrix H(T, R);

        return kin->getInverseConfiguration(std::move(currentJoints), H);
    }
}

void dxKinInterface::test()
{
    std::vector<double> q(7, 0);

    q[1] = -90;
    q[3] = -50;

    vpHomogeneousMatrix H = this->kin->getForwardPosition(q);

    vpTranslationVector T(H);
    vpRotationMatrix R(H);
    vpRzyxVector rv(R);


    //-------------------------------
    //		Inverse kinematics from Kin
    //-------------------------------

    vpTranslationVector T_tar;
    T_tar[0] = T[0] + 0.15;
    T_tar[1] = T[1] + 0.30;
    T_tar[2] = T[2] + 0.20;

    vpRotationMatrix R_tar;
    R_tar.buildFrom(H);
    vpRzyxVector rv_tar(R_tar);

    vpHomogeneousMatrix H_tar(T_tar, R_tar);


    auto start = std::chrono::high_resolution_clock::now();

    std::vector<double> r1 = this->kin->getInverseConfiguration(q, H_tar);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);

    std::cout << std::fixed;
    std::cout << "Seed: ";
    printVector(q);

    std::cout << "R1: ";
    printVector(r1);

    vpHomogeneousMatrix H1 = this->kin->getForwardPosition(r1);
    vpTranslationVector T1 = H1.getTranslationVector();
    vpRzyxVector rv1(H1.getRotationMatrix());

    std::cout << std::fixed;
    std::cout << "Target: " << T_tar[0] << "  " << T_tar[1] << "  " << T_tar[2] << "  " << vpMath::deg(rv_tar[0]) << "  " << vpMath::deg(rv_tar[1]) << "  " << vpMath::deg(rv_tar[2]) << std::endl;
    std::cout << "  IK : " << T1[0] << "  " << T1[1] << "  " << T1[2] << "  " << vpMath::deg(rv1[0]) << "  " << vpMath::deg(rv1[1]) << "  " << vpMath::deg(rv1[2]) << std::endl;

    std::cout << "  Duration: " << duration.count() << " \xE6s " << std::endl;



    //-------------------------------
    //		Inverse kinematics within Function
    //-------------------------------
    std::cout << "  -------------------------------------------- " << std::endl;

    vpQuaternionVector Q(R_tar);

    std::vector<double> tp(7, 0);
    tp[0] = T_tar[0];
    tp[1] = T_tar[1];
    tp[2] = T_tar[2];
    tp[3] = Q.x();
    tp[4] = Q.y();
    tp[5] = Q.z();
    tp[6] = Q.w();



    auto strt = std::chrono::high_resolution_clock::now();

    std::vector<double> r2 = this->getJointCommand(q, tp);

    auto stp = std::chrono::high_resolution_clock::now();
    auto dur = std::chrono::duration_cast <std::chrono::microseconds> (stp - strt);

    std::cout << std::fixed;
    std::cout << "q: ";
    printVector(q);

    std::cout << "R2: ";
    printVector(r2);

    vpHomogeneousMatrix H2 = this->kin->getForwardPosition(r2);
    vpTranslationVector T2 = H2.getTranslationVector();
    vpRzyxVector rv2(H2.getRotationMatrix());

    std::cout << std::fixed;
    std::cout << "Target: " << T_tar[0] << "  " << T_tar[1] << "  " << T_tar[2] << "  " << vpMath::deg(rv_tar[0]) << "  " << vpMath::deg(rv_tar[1]) << "  " << vpMath::deg(rv_tar[2]) << std::endl;
    std::cout << "  IK : " << T2[0] << "  " << T2[1] << "  " << T2[2] << "  " << vpMath::deg(rv2[0]) << "  " << vpMath::deg(rv2[1]) << "  " << vpMath::deg(rv2[2]) << std::endl;

    std::cout << "  Duration: " << dur.count() << " \xE6s " << std::endl;
}
