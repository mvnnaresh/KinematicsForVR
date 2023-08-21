// kinematicsForVR.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <chrono>
#include "dxKin.h"

using namespace std;


template <typename T>
void printVector(vector<T> v)
{
    for_each(v.begin(), v.end(), [](T x)
    {
        cout << x << "   ";
    });
    cout << endl;
}


int main()
{
    dxKin kin;

    vector<double> q(kin.numJoints, 0);

    q[1] = -90;
    q[3] = -50;

    vpHomogeneousMatrix H = kin.getForwardPosition(q);

    vpTranslationVector T(H);
    vpRotationMatrix R(H);
    vpRzyxVector rv(R);


    //-------------------------------
    //		Inverse kinematics
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

    vector<double> r1 = kin.getJointCommand(q, H_tar);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);

    cout << std::fixed;
    cout << "Seed: ";
    printVector(q);

    cout << "R1: ";
    printVector(r1);

    vpHomogeneousMatrix H1 = kin.getForwardPosition(r1);
    vpTranslationVector T1 = H1.getTranslationVector();
    vpRzyxVector rv1(H1.getRotationMatrix());

    cout << std::fixed;
    cout << "Target: " << T_tar[0] << "  " << T_tar[1] << "  " << T_tar[2] << "  " << vpMath::deg(rv_tar[0]) << "  " << vpMath::deg(rv_tar[1]) << "  " << vpMath::deg(rv_tar[2]) << endl;
    cout << "  IK : " << T1[0] << "  " << T1[1] << "  " << T1[2] << "  " << vpMath::deg(rv1[0]) << "  " << vpMath::deg(rv1[1]) << "  " << vpMath::deg(rv1[2]) << endl;

    cout << "  KDL: " << duration.count() << " \xE6s " << endl;
}

