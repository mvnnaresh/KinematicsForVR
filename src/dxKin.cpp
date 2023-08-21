#include "dxKin.h"


dxKin::dxKin(ChainData robotData) :
    IK_METHOD(static_cast<IKMethod>(100)),
    fksolver(nullptr),
    iksolver1v(nullptr),
    iksolver(nullptr)
{
    this->robot = robotData;

    this->low = this->robot.q_min;
    this->high = this->robot.q_max;
    this->velLim = this->robot.q_vmax;
    this->numJoints = this->robot.chain.getNrOfJoints();

    //---- Safety offset for joints to not reach close to limits
    double SAFTETY_OFFSET = 1 * M_PI / 180;

    //---- Allocate the memory for the joints position, velocity and acceleration
    q = std::make_shared<KDL::JntArray>(this->numJoints);
    qdot = std::make_shared<KDL::JntArray>(this->numJoints);
    qdotdot = std::make_shared<KDL::JntArray>(this->numJoints);

    //---- Creation of joint arrays:
    for (auto i = 0; i < numJoints; i++)
    {
        (*q)(i) = 0.0;
        (*qdot)(i) = 0.0;
        (*qdotdot)(i) = 0.0;
    }

    //---- load joint limits to KDL joint arrays
    KDL::JntArray q_min(this->numJoints);
    KDL::JntArray q_max(this->numJoints);
    for (int i = 0; i < this->numJoints; i++)
    {
        q_min(i) = this->robot.q_min[i] + SAFTETY_OFFSET;
        q_max(i) = this->robot.q_max[i] + SAFTETY_OFFSET;
    }

    //---- Kinematic solvers of KDL
    fksolver = new KDL::ChainFkSolverPos_recursive(this->robot.chain);
    iksolver1v = new KDL::ChainIkSolverVel_pinv(this->robot.chain);
    iksolver = new KDL::ChainIkSolverPos_NR_JL(this->robot.chain, q_min, q_max, *fksolver, *iksolver1v, 1000, 1e-4);

}

std::vector<double> dxKin::getInverseConfiguration(std::vector<double> qcurr, vpHomogeneousMatrix H)
{
    if (IK_METHOD == IKMethod::KDL_PINV_IK)
    {
        return this->ikSolverKDL(H, qcurr);
    }
    if (IK_METHOD == IKMethod::JPINV_NSO)
    {
        return this->ikSolverNSO(H, qcurr);
    }
    if (IK_METHOD == IKMethod::DLS_NSO)
    {
        return this->ikSolverDLSNSO(H, qcurr);
    }
    else
    {
        std::cerr << " Invalid IK method. Returning seed!!" << std::endl;
        return qcurr;
    }
}

void dxKin::setIKMethod(IKMethod method)
{
    this->IK_METHOD = method;
}

//------------------------------------------------------------
//					KINEMATIC MEMBERS
//------------------------------------------------------------
vpHomogeneousMatrix dxKin::getForwardPosition(std::vector <double> j)
{
    KDL::Frame F = this->fkSolver(j);

    return KDLFrameTovpHomMat(F);
}

void dxKin::setJointPos(std::vector<double> pos)
{
    assert(pos.size() == numJoints);

    for (unsigned int i = 0; i < this->numJoints; i++)
    {
        (*q)(i) = RADIANS(pos[i]);
    }
}

//------------------------------------------------------------
//					JACOBIAN MATRIX
//------------------------------------------------------------
vpMatrix dxKin::getBaseJacobianMatrix(std::vector<double> q)
{
    vpMatrix J_vp;

    KDL::JntArray jointAngles(this->numJoints);
    // Assign joint positions
    for (unsigned int i = 0; i < this->numJoints; i++)
    {
        jointAngles(i) = (double)vpMath::rad(q[i]);
    }

    KDL::ChainJntToJacSolver jacobianSolver(this->robot.chain);
    KDL::Jacobian baseJacobian(this->numJoints);

    // Calculate the Jacobian at the base frame given the joint angles
    int result = jacobianSolver.JntToJac(jointAngles, baseJacobian);

    if (result < 0)
    {
        std::cerr << "Failed to calculate the Jacobian!\n";
    }

    J_vp.resize(baseJacobian.rows(), baseJacobian.columns(), true);

    // Assign Jacobian matrix
    for (unsigned int i = 0; i < baseJacobian.rows(); ++i)
    {
        for (unsigned int j = 0; j < baseJacobian.columns(); ++j)
        {
            J_vp[i][j] = baseJacobian(i, j);
        }
    }

    return J_vp;
}

//------------------------------------------------------------
//					KINEMATIC SOLVERS
// Note: All INPUT and OUTPUT angles are in DEGREES
//------------------------------------------------------------

KDL::Frame dxKin::fkSolver(std::vector <double> j)
{
    this->setJointPos(j);

    KDL::Frame cartPos;
    bool kinematics_status = this->fksolver->JntToCart(*this->q, cartPos);

    return cartPos;
}

std::vector <double> dxKin::ikSolverKDL(vpHomogeneousMatrix cartPos, std::vector<double> seed)
{
    std::vector<double> q_sol;
    q_sol.assign(this->numJoints, 0);

    //assign joint values - create KDL arrays
    KDL::JntArray q_init(this->numJoints);
    KDL::JntArray q_out(this->numJoints);

    //assign seed values
    for (int i = 0; i < this->numJoints; i++)
        q_init(i) = RADIANS(seed[i]);

    vpRotationMatrix R(cartPos);
    vpTranslationVector T(cartPos);
    KDL::Frame F_in;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
            F_in.M(i, j) = R[i][j];
    }

    for (int i = 0; i < 3; i++)
        F_in.p(i) = T[i];

    int ret = this->iksolver->CartToJnt(q_init, F_in, q_out);

    // Interprete results
    if (ret < 0) //Some problem in IK solution
    {
        std::cout << "KDL IK solver FAILED - No solution found. Returning seed configuration!!" << std::endl;
        q_sol = seed;
    }
    else
    {
        //dxUtils::printMessage2("IK solver SUCCESS - solution found!!");
        for (int i = 0; i < numJoints; i++)
            q_sol[i] = q_out(i)*(180 / M_PI);
    }
    return q_sol;
}

std::vector <double> dxKin::ikSolverNSO(vpHomogeneousMatrix cartPos, std::vector<double> seed)
{
    vpMatrix JPInv;
    vpColVector q_curr = seed;

    int numDOF = this->numJoints;
    int iterations = 0;
    int maxIterations = 1000;		// Maximum number of iterations to avoid infinite loops

    double tolerance = 1e-5;		// Tolerance for convergence
    double manipulabilityPrev = 0;
    double manipulabilityGain = 1.2;

    std::vector<double> q_prev(numDOF, 0);
    std::vector<double> NSq0_des(numDOF, 0);
    std::vector<double> NSq0_dot(numDOF, 0);

    KDL::Frame targetFrame = vpHomMatToKDLFrame(cartPos);

    while (iterations < maxIterations)
    {
        vpHomogeneousMatrix H = this->getForwardPosition(q_curr.toStdVector());

        KDL::Frame currentFrame = vpHomMatToKDLFrame(H);

        if (KDL::Equal(currentFrame, targetFrame, tolerance))
        {
            break;
        }

        vpColVector positionError = this->computePoseError(cartPos, H);

        //Nullspace control
        double manipulability = this->getManipulability(q_curr.toStdVector());

        if (iterations > 0)
        {
            // Applying secondary objective for maximising the manipulability Eq.3.56
            for (int i = 0; i < numDOF; i++)
            {
                //Computing the delta of join pose for Eq.3.55
                double dq = (RADIANS(q_curr[i] - q_prev[i]));
                if (fabs(dq) < 1e-6)
                {
                    NSq0_dot[i] = 0;
                    continue;
                }
                //Needs to be multiplied by dt because of the integration
                NSq0_dot[i] = manipulabilityGain * (manipulability - manipulabilityPrev) / dq;
            }
        }

        //Get the nullspace command
        std::vector<double> qcmd = this->getPositionWithNullSpaceAndDLS(positionError.toStdVector(), q_curr.toStdVector(), NSq0_dot);

        q_curr += vpColVector(qcmd);

        q_prev = q_curr.toStdVector();
        manipulabilityPrev = manipulability;

        iterations++;
    }
    if (iterations >= maxIterations)
    {
        std::cout << "KDL IK solver FAILED - No solution found. Returning seed configuration!!" << std::endl;
    }

    return q_curr.toStdVector();
}

std::vector <double> dxKin::ikSolverDLSNSO(vpHomogeneousMatrix cartPos, std::vector<double> seed)
{
    vpMatrix JPInv;
    vpColVector q_curr = seed;

    int numDOF = this->numJoints;
    int iterations = 0;
    int maxIterations = 1000;		// Maximum number of iterations to avoid infinite loops

    double tolerance = 1e-5;		// Tolerance for convergence
    double manipulabilityPrev = 0;
    double manipulabilityGain = 1.2;

    std::vector<double> q_prev(numDOF, 0);
    std::vector<double> NSq0_des(numDOF, 0);
    std::vector<double> NSq0_dot(numDOF, 0);

    KDL::Frame targetFrame = vpHomMatToKDLFrame(cartPos);

    while (iterations < maxIterations)
    {
        vpHomogeneousMatrix H = this->getForwardPosition(q_curr.toStdVector());

        KDL::Frame currentFrame = vpHomMatToKDLFrame(H);

        if (KDL::Equal(currentFrame, targetFrame, tolerance))
        {
            break;
        }

        vpColVector positionError = this->computePoseError(cartPos, H);

        //Nullspace control
        double manipulability = this->getManipulability(q_curr.toStdVector());

        if (iterations > 0)
        {
            // Applying secondary objective for maximising the manipulability Eq.3.56
            for (int i = 0; i < numDOF; i++)
            {
                //Computing the delta of join pose for Eq.3.55
                double dq = (RADIANS(q_curr[i] - q_prev[i]));
                if (fabs(dq) < 1e-6)
                {
                    NSq0_dot[i] = 0;
                    continue;
                }
                //Needs to be multiplied by dt because of the integration
                NSq0_dot[i] = manipulabilityGain * (manipulability - manipulabilityPrev) / dq;
            }
        }

        //Get the nullspace command
        std::vector<double> qcmd = this->getPositionWithNullSpace(positionError.toStdVector(), q_curr.toStdVector(), NSq0_dot);

        q_curr += vpColVector(qcmd);

        q_prev = q_curr.toStdVector();
        manipulabilityPrev = manipulability;

        iterations++;
    }
    if (iterations >= maxIterations)
    {
        std::cout << "KDL IK solver FAILED - No solution found. Returning seed configuration!!" << std::endl;
    }

    return q_curr.toStdVector();
}
//------------------------------------------------------------
//			CONTROL BY NULLSPACE OPTIMISATION
//------------------------------------------------------------
std::vector<double> dxKin::getPositionWithNullSpace(std::vector<double> xdot, std::vector<double> qcurrent, std::vector<double> nullspaceQ0_dot)
{
    vpMatrix Jpinv;
    int numDOF = this->numJoints;
    vpMatrix I;
    I.eye(numDOF, numDOF);

    vpColVector q_curr = qcurrent;
    vpColVector q0_dot = nullspaceQ0_dot;
    vpColVector x_dot = xdot;

    //Jacobian matrix
    vpMatrix J = this->getBaseJacobianMatrix(q_curr.toStdVector());

    //Jacobian Pseudo Inverse
    Jpinv = J.pseudoInverse();

    //q0_dot is the secondary task that needs to be optimised
    vpMatrix Nspace = (I - (Jpinv * J))*q0_dot;

    //Sciciliano P.125 - Eq. 3.54
    //WARNING: This equation only works when the Jacobian is fullrank!!
    vpColVector q_cmd_dot = Jpinv * x_dot + Nspace;

    return q_cmd_dot.toStdVector();
}

std::vector<double> dxKin::getPositionWithNullSpaceAndDLS(std::vector<double> xdot, std::vector<double> qcurrent, std::vector<double> nullspaceQ0_dot)
{
    vpMatrix Jpinv;
    int numDOF = this->numJoints;
    vpMatrix I;
    I.eye(numDOF, numDOF);

    vpMatrix I_LM;
    I_LM.eye(6);
    double lambda = DLS_LAMBDA;

    vpColVector q_curr = qcurrent;
    vpColVector q0_dot = nullspaceQ0_dot;
    vpColVector x_dot = xdot;

    //Jacobian matrix
    vpMatrix J = this->getBaseJacobianMatrix(q_curr.toStdVector());

    //Jacobian with damped least squares
    vpMatrix J_LM = J.transpose() * ((J * J.transpose()) + (std::pow(lambda, 2) * I_LM)).inverseByLU();

    //q0_dot is the secondary task that needs to be optimised
    vpMatrix Nspace = (I - (J_LM * J))*q0_dot;

    //Sciciliano P.125 - Eq. 3.54
    //WARNING: This equation only works when the Jacobian is fullrank!!
    vpColVector q_cmd = J_LM * x_dot + Nspace + q_curr;

    return q_cmd.toStdVector();
}

//------------------------------------------------------------
//				ERROR COMPUTATION FUNCTIONS
//------------------------------------------------------------
vpColVector dxKin::computePoseError(vpHomogeneousMatrix Hd, vpHomogeneousMatrix Hc)
{
    vpColVector ret(6, 0);

    vpTranslationVector T1 = Hd.getTranslationVector();
    vpTranslationVector T2 = Hc.getTranslationVector();

    vpQuaternionVector Q1, Q2;
    Q1.buildFrom(Hd.getRotationMatrix());
    Q2.buildFrom(Hc.getRotationMatrix());

    vpColVector qe = this->computeQuaternionError(Q1, Q2);

    ret[0] = T1[0] - T2[0];
    ret[1] = T1[1] - T2[1];
    ret[2] = T1[2] - T2[2];
    ret[3] = qe[0];
    ret[4] = qe[1];
    ret[5] = qe[2];

    return ret;
}

vpColVector dxKin::computeQuaternionError(vpQuaternionVector Qd, vpQuaternionVector Qc)
{
    double sd = Qd.w(), sc = Qc.w();

    vpColVector ud(3), uc(3);
    ud[0] = Qd.x();
    ud[1] = Qd.y();
    ud[2] = Qd.z();

    uc[0] = Qc.x();
    uc[1] = Qc.y();
    uc[2] = Qc.z();

    vpColVector e(3);
    //e = sd*uc - sc*ud + vpColVector::skew(ud)*uc;
    //e = -sd*uc + sc*ud - vpColVector::skew(ud)*uc;
    double e_s = sc * sd + vpColVector::dotProd(ud, uc);

    e = sc * ud - sd * uc - vpColVector::cross(ud, uc);

    e = round(e_s) * e;

    return e;
}

double dxKin::getManipulability(std::vector<double> qcurrent)
{
    vpMatrix I;
    I.eye(this->numJoints, this->numJoints);

    vpMatrix J = this->getBaseJacobianMatrix(qcurrent);

    return std::sqrt((J * J.transpose()).det());
}
