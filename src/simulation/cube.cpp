#include "cube.h"

#include "Eigen/Dense"

#include "../util/helper.h"
namespace simulation {
constexpr float g_cdK = 2000.0f;
constexpr float g_cdD = 50.0f;

Cube::Cube()
    : particleNumPerEdge(10),
      cubeLength(2.0),
      initialPosition(Eigen::Vector3f(0.0, 0.0, 0.0)),
      springCoefStruct(g_cdK),
      springCoefShear(g_cdK),
      springCoefBending(g_cdK),
      damperCoefStruct(g_cdD),
      damperCoefShear(g_cdD),
      damperCoefBending(g_cdD) {
    particleNumPerFace = particleNumPerEdge * particleNumPerEdge;
    initializeParticle();
    initializeSpring();
}

Cube::Cube(const Eigen::Vector3f &a_kInitPos, const float cubeLength, const int numAtEdge, const float dSpringCoef,
           const float dDamperCoef)
    : particleNumPerEdge(numAtEdge),
      cubeLength(cubeLength),
      initialPosition(a_kInitPos),
      springCoefStruct(dSpringCoef),
      springCoefShear(dSpringCoef),
      springCoefBending(dSpringCoef),
      damperCoefStruct(dDamperCoef),
      damperCoefShear(dDamperCoef),
      damperCoefBending(dDamperCoef) {
    particleNumPerFace = numAtEdge * numAtEdge;
    initializeParticle();
    initializeSpring();
}

int Cube::getParticleNum() const { return static_cast<int>(particles.size()); }

int Cube::getSpringNum() const { return static_cast<int>(springs.size()); }

int Cube::getNumAtEdge() const { return particleNumPerEdge; }

unsigned int Cube::getPointMap(const int a_ciSide, const int a_ciI, const int a_ciJ) {
    int r = -1;

    switch (a_ciSide) {
        case 1:  // [a_ciI][a_ciJ][0] bottom face
            r = particleNumPerFace * a_ciI + particleNumPerEdge * a_ciJ;
            break;
        case 6:  // [a_ciI][a_ciJ][9] top face
            r = particleNumPerFace * a_ciI + particleNumPerEdge * a_ciJ + particleNumPerEdge - 1;
            break;
        case 2:  // [a_ciI][0][a_ciJ] front face
            r = particleNumPerFace * a_ciI + a_ciJ;
            break;
        case 5:  // [a_ciI][9][a_ciJ] back face
            r = particleNumPerFace * a_ciI + particleNumPerEdge * (particleNumPerEdge - 1) + a_ciJ;
            break;
        case 3:  // [0][a_ciI][a_ciJ] left face
            r = particleNumPerEdge * a_ciI + a_ciJ;
            break;
        case 4:  // [9][a_ciI][a_ciJ] ra_ciIght face
            r = particleNumPerFace * (particleNumPerEdge - 1) + particleNumPerEdge * a_ciI + a_ciJ;
            break;
    }

    return r;
}

Particle &Cube::getParticle(int particleIdx) { return particles[particleIdx]; }

std::vector<Particle> *Cube::getParticlePointer() { return &particles; }

Spring &Cube::getSpring(int springIdx) { return springs[springIdx]; }

void Cube::setSpringCoef(const float springCoef, const Spring::SpringType springType) {
    if (springType == Spring::SpringType::STRUCT) {
        springCoefStruct = springCoef;
        updateSpringCoef(springCoef, Spring::SpringType::STRUCT);
    } else if (springType == Spring::SpringType::SHEAR) {
        springCoefShear = springCoef;
        updateSpringCoef(springCoef, Spring::SpringType::SHEAR);
    } else if (springType == Spring::SpringType::BENDING) {
        springCoefBending = springCoef;
        updateSpringCoef(springCoef, Spring::SpringType::BENDING);
    }
}

void Cube::setDamperCoef(const float damperCoef, const Spring::SpringType springType) {
    if (springType == Spring::SpringType::STRUCT) {
        damperCoefStruct = damperCoef;
        updateDamperCoef(damperCoef, Spring::SpringType::STRUCT);
    } else if (springType == Spring::SpringType::SHEAR) {
        damperCoefShear = damperCoef;
        updateDamperCoef(damperCoef, Spring::SpringType::SHEAR);
    } else if (springType == Spring::SpringType::BENDING) {
        damperCoefBending = damperCoef;
        updateDamperCoef(damperCoef, Spring::SpringType::BENDING);
    }
}

void Cube::resetCube(const Eigen::Vector3f &offset, const float &rotate) {
    float dTheta = util::radians(rotate);  //  change angle from degree to
                                           //  radian

    for (unsigned int uiI = 0; uiI < particles.size(); uiI++) {
        int i = uiI / particleNumPerFace;
        int j = (uiI / particleNumPerEdge) % particleNumPerEdge;
        int k = uiI % particleNumPerEdge;
        float offset_x = (float)((i - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));
        float offset_y = (float)((j - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));
        float offset_z = (float)((k - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));

        Eigen::Vector3f RotateVec(offset_x, offset_y,
                                  offset_z);  //  vector from center of cube to the particle

        Eigen::AngleAxis<float> rotation(dTheta, Eigen::Vector3f(1.0f, 0.0f, 1.0f).normalized());

        RotateVec = rotation * RotateVec;

        particles[uiI].setPosition(initialPosition + offset + RotateVec);
        particles[uiI].setForce(Eigen::Vector3f::Zero());
        particles[uiI].setVelocity(Eigen::Vector3f::Zero());
    }
}

void Cube::addForceField(const Eigen::Vector3f &force) {
    for (unsigned int uiI = 0; uiI < particles.size(); uiI++) {
        particles[uiI].setAcceleration(force);
    }
}

void Cube::computeInternalForce() {

    // TODO
    // Trace every spring and apply the force accordingly

    for (int i = 0; i < springs.size(); i++) {

        // Compute spring force and damper force for every spring
        // Add forces to both particles on the two ends of the spring

        Spring Spring = springs[i];

        Particle StartParticle = particles[Spring.getSpringStartID()];
        Particle EndParticle = particles[Spring.getSpringEndID()];

        Eigen::Vector3f SpringForce = computeSpringForce(StartParticle.getPosition(), EndParticle.getPosition(),
                                                         Spring.getSpringCoef(), Spring.getSpringRestLength());

        Eigen::Vector3f DamperForce = computeDamperForce(StartParticle.getPosition(), EndParticle.getPosition(),
                                                         StartParticle.getVelocity(), EndParticle.getVelocity(),
                                                         Spring.getDamperCoef());

        StartParticle.addForce(SpringForce);
        StartParticle.addForce(DamperForce);
        particles[Spring.getSpringStartID()] = StartParticle;

        EndParticle.addForce(-SpringForce);
        EndParticle.addForce(-DamperForce);
        particles[Spring.getSpringEndID()] = EndParticle;
    }
}

Eigen::Vector3f Cube::computeSpringForce(const Eigen::Vector3f &positionA, const Eigen::Vector3f &positionB,
                                         const float springCoef, const float restLength) {
    
    float epsilon = 0.001;
    
    // TODO
    // f_s = -k_s(|x_a - x_b| - r) * l, k_s > 0

    Eigen::Vector3f relativePosition = positionA - positionB;
    float currentLength = relativePosition.norm();

    Eigen::Vector3f normalizedRelativePosition = relativePosition / currentLength;

    Eigen::Vector3f SpringForce = -springCoef * ((currentLength - restLength)) * normalizedRelativePosition;

    return SpringForce;
}

Eigen::Vector3f Cube::computeDamperForce(const Eigen::Vector3f &positionA, const Eigen::Vector3f &positionB,
                                         const Eigen::Vector3f &velocityA, const Eigen::Vector3f &velocityB,
                                         const float damperCoef) {
    // TODO
    // f_d = -k_d((v_a - v_b) . l) * l, k_d > 0

    Eigen::Vector3f relativePosition = positionA - positionB;
    Eigen::Vector3f relativeVelocity = velocityA - velocityB;
    float currentLength = relativePosition.norm();

    Eigen::Vector3f normalizedRelativePosition = relativePosition / currentLength;

    Eigen::Vector3f DamperForce = -damperCoef * (relativeVelocity.dot(normalizedRelativePosition)) * normalizedRelativePosition;

    return DamperForce;
}

void Cube::initializeParticle() {
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                Particle Particle;
                float offset_x = (float)((i - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));
                float offset_y = (float)((j - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));
                float offset_z = (float)((k - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));
                Particle.setPosition(Eigen::Vector3f(initialPosition(0) + offset_x, initialPosition(1) + offset_y,
                                                     initialPosition(2) + offset_z));
                particles.push_back(Particle);
            }
        }
    }
}

void Cube::initializeSpring() {

    int iParticleID = 0;
    int iNeighborID = 0;
    Eigen::Vector3f SpringStartPos;
    Eigen::Vector3f SpringEndPos;
    float Length;

    float springCoef = g_cdK;
    float damperCoef = g_cdD;

    // TODO
    // Initialize three types of springs
    // Structure, bending, and shear springs

    // Structure springs connect the center particle to the 6 particles on 3 axes (2 directions each).
    // Bending springs connect the center particle to the 6 particles on 3 axes (2 directions each) 
    //     that are 2 particles away from it to avoid structure springs from bending.
    // Shear springs connect the center particle to all of its diagonal particles. (20 of them)
    //     (26 neighboring pariticles - 6 pariticles on the axes)

    // Each spring connects 2 particles ===> Each particle handles half of the springs

    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {

                iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;

                // Struct Springs 
                // Each particle is responsible for only the positive directions
                // i.e. at most 3 struct springs for each particle

                std::vector<int> Struct_x;
                std::vector<int> Struct_y;
                std::vector<int> Struct_z;

                if (i < particleNumPerEdge - 1) {
                    Struct_x.push_back(i + 1);
                    Struct_y.push_back(j);
                    Struct_z.push_back(k);
                }

                if (j < particleNumPerEdge - 1) {
                    Struct_x.push_back(i);
                    Struct_y.push_back(j + 1);
                    Struct_z.push_back(k);
                }

                if (k < particleNumPerEdge - 1) {
                    Struct_x.push_back(i);
                    Struct_y.push_back(j);
                    Struct_z.push_back(k + 1);
                }

                for (int t = 0; t < Struct_x.size(); t++) {
                    iNeighborID = Struct_x[t] * particleNumPerFace + Struct_y[t] * particleNumPerEdge + Struct_z[t];

                    SpringStartPos = particles[iParticleID].getPosition();
                    SpringEndPos = particles[iNeighborID].getPosition();
                    Length = (SpringStartPos - SpringEndPos).norm();

                    Spring Spring(iParticleID, iNeighborID, Length, springCoef, damperCoef, Spring::SpringType::STRUCT);
                    springs.push_back(Spring);
                }

                // Bend Springs
                // Each particle is responsible for only the positive directions
                // i.e. at most 3 bend springs for each particle

                std::vector<int> Bend_x;
                std::vector<int> Bend_y;
                std::vector<int> Bend_z;

                if (i < particleNumPerEdge - 2) {
                    Bend_x.push_back(i + 2);
                    Bend_y.push_back(j);
                    Bend_z.push_back(k);
                }

                if (j < particleNumPerEdge - 2) {
                    Bend_x.push_back(i);
                    Bend_y.push_back(j + 2);
                    Bend_z.push_back(k);
                }

                if (k < particleNumPerEdge - 2) {
                    Bend_x.push_back(i);
                    Bend_y.push_back(j);
                    Bend_z.push_back(k + 2);
                }
                
                if (i == 0) {
                    Bend_x.push_back(particleNumPerEdge - 1);
                    Bend_y.push_back(j);
                    Bend_z.push_back(k);
                }

                if (j == 0) {
                    Bend_x.push_back(i);
                    Bend_y.push_back(particleNumPerEdge - 1);
                    Bend_z.push_back(k);
                }

                if (k == 0) {
                    Bend_x.push_back(i);
                    Bend_y.push_back(j);
                    Bend_z.push_back(particleNumPerEdge - 1);
                }
                
                for (int t = 0; t < Bend_x.size(); t++) {
                    iNeighborID = Bend_x[t] * particleNumPerFace + Bend_y[t] * particleNumPerEdge + Bend_z[t];

                    SpringStartPos = particles[iParticleID].getPosition();
                    SpringEndPos = particles[iNeighborID].getPosition();
                    Length = (SpringStartPos - SpringEndPos).norm();

                    Spring Spring(iParticleID, iNeighborID, Length, springCoef, damperCoef, Spring::SpringType::BENDING);
                    springs.push_back(Spring);
                }

                // Shear Springs
                // Each particle is responsible for 10 shear springs

                std::vector<int> Shear_x;
                std::vector<int> Shear_y;
                std::vector<int> Shear_z;

                // Directions with all positive axes

                if (i < particleNumPerEdge - 1 && j < particleNumPerEdge - 1) {
                    Shear_x.push_back(i + 1);
                    Shear_y.push_back(j + 1);
                    Shear_z.push_back(k);
                }

                if (j < particleNumPerEdge - 1 && k < particleNumPerEdge - 1) {
                    Shear_x.push_back(i);
                    Shear_y.push_back(j + 1);
                    Shear_z.push_back(k + 1);
                }

                if (i < particleNumPerEdge - 1 && k < particleNumPerEdge - 1) {
                    Shear_x.push_back(i + 1);
                    Shear_y.push_back(j);
                    Shear_z.push_back(k + 1);
                }

                // Directions with negative axis
                
                if (i > 0 && j < particleNumPerEdge - 1) {
                    Shear_x.push_back(i - 1);
                    Shear_y.push_back(j + 1);
                    Shear_z.push_back(k);
                }

                if (j > 0 && k < particleNumPerEdge - 1) {
                    Shear_x.push_back(i);
                    Shear_y.push_back(j - 1);
                    Shear_z.push_back(k + 1);
                }

                if (i > 0 && k < particleNumPerEdge - 1) {
                    Shear_x.push_back(i - 1);
                    Shear_y.push_back(j);
                    Shear_z.push_back(k + 1);
                }

                // Directions that are cubical diagonal

                if (i < particleNumPerEdge - 1 && j < particleNumPerEdge - 1 && k < particleNumPerEdge - 1) {
                    Shear_x.push_back(i + 1);
                    Shear_y.push_back(j + 1);
                    Shear_z.push_back(k + 1);
                }

                if (i > 0 && j < particleNumPerEdge - 1 && k < particleNumPerEdge - 1) {
                    Shear_x.push_back(i - 1);
                    Shear_y.push_back(j + 1);
                    Shear_z.push_back(k + 1);
                }

                if (i < particleNumPerEdge - 1 && j > 0 && k < particleNumPerEdge - 1) {
                    Shear_x.push_back(i + 1);
                    Shear_y.push_back(j - 1);
                    Shear_z.push_back(k + 1);
                }

                if (i > 0 && j > 0 && k < particleNumPerEdge - 1) {
                    Shear_x.push_back(i - 1);
                    Shear_y.push_back(j - 1);
                    Shear_z.push_back(k + 1);
                }

                for (int t = 0; t < Shear_x.size(); t++) {
                    iNeighborID = Shear_x[t] * particleNumPerFace + Shear_y[t] * particleNumPerEdge + Shear_z[t];

                    SpringStartPos = particles[iParticleID].getPosition();
                    SpringEndPos = particles[iNeighborID].getPosition();
                    Length = (SpringStartPos - SpringEndPos).norm();

                    Spring Spring(iParticleID, iNeighborID, Length, springCoef, damperCoef, Spring::SpringType::SHEAR);
                    springs.push_back(Spring);
                }
            }
        }
    }
}

void Cube::updateSpringCoef(const float a_cdSpringCoef, const Spring::SpringType a_cSpringType) {
    for (unsigned int uiI = 0; uiI < springs.size(); uiI++) {
        if (springs[uiI].getType() == a_cSpringType) {
            springs[uiI].setSpringCoef(a_cdSpringCoef);
        }
    }
}

void Cube::updateDamperCoef(const float a_cdDamperCoef, const Spring::SpringType a_cSpringType) {
    for (unsigned int uiI = 0; uiI < springs.size(); uiI++) {
        if (springs[uiI].getType() == a_cSpringType) {
            springs[uiI].setDamperCoef(a_cdDamperCoef);
        }
    }
}
}  //  namespace simulation
