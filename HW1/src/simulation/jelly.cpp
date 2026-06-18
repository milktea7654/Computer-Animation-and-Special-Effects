#include "jelly.h"

#include "Eigen/Dense"

#include "../util/helper.h"
#include <iostream>
namespace simulation {
constexpr float g_cdK = 2500.0f;
constexpr float g_cdD = 50.0f;

Jelly::Jelly()
    : particleNumPerEdge(10),
      jellyLength(2.0),
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

Jelly::Jelly(const Eigen::Vector3f &a_kInitPos, const float jellyLength, const int numAtEdge, const float dSpringCoef,
           const float dDamperCoef)
    : particleNumPerEdge(numAtEdge),
      jellyLength(jellyLength),
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

int Jelly::getParticleNum() const { return static_cast<int>(particles.size()); }

int Jelly::getSpringNum() const { return static_cast<int>(springs.size()); }

int Jelly::getNumAtEdge() const { return particleNumPerEdge; }

unsigned int Jelly::getPointMap(const int a_ciSide, const int a_ciI, const int a_ciJ) {
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

Particle &Jelly::getParticle(int particleIdx) { return particles[particleIdx]; }

std::vector<Particle> *Jelly::getParticlePointer() { return &particles; }

Spring &Jelly::getSpring(int springIdx) { return springs[springIdx]; }

void Jelly::setSpringCoef(const float springCoef, const Spring::SpringType springType) {
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

void Jelly::setDamperCoef(const float damperCoef, const Spring::SpringType springType) {
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

void Jelly::resetJelly(const Eigen::Vector3f &offset, const float &rotate) {
    float dTheta = util::radians(rotate);  //  change angle from degree to
                                           //  radian

    for (unsigned int uiI = 0; uiI < particles.size(); uiI++) {
        int i = uiI / particleNumPerFace;
        int j = (uiI / particleNumPerEdge) % particleNumPerEdge;
        int k = uiI % particleNumPerEdge;
        float offset_x = (float)((i - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));
        float offset_y = (float)((j - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));
        float offset_z = (float)((k - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));

        Eigen::Vector3f RotateVec(offset_x, offset_y,
                                  offset_z);  //  vector from center of cube to the particle

        Eigen::AngleAxis<float> rotation(dTheta, Eigen::Vector3f(1.0f, 0.0f, 1.0f).normalized());

        RotateVec = rotation * RotateVec;

        particles[uiI].setPosition(initialPosition + offset + RotateVec);
        particles[uiI].setForce(Eigen::Vector3f::Zero());
        particles[uiI].setVelocity(Eigen::Vector3f::Zero());
    }
}

void Jelly::addForceField(const Eigen::Vector3f &force) {
    for (unsigned int uiI = 0; uiI < particles.size(); uiI++) {
        particles[uiI].setAcceleration(force);
    }
}

void Jelly::computeInternalForce()
{
    // TODO#2-3: Compute the internal force (including spring force and damper force) for each spring.
    //   1. Read the start-particle and end-particle index from spring.
    //   2. Use `getPosition()` to get particle i's position and `getVelocity()` to get particle i's velocity.
    //   3. Call `computeSpringForce` and `computeDamperForce` to compute spring force and damper force.
    //   4. Compute net internal force and call `addForce` to apply the force onto particles.
    // Note:
    //   1. Direction of the force.
    // Hint:
    //   1. Use a.norm() to get length of a
    //   2. Use a.normalize() to normalize a inplace.
    //      a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.

    for (int i = 0; i < springs.size(); i++) {
        int springStartID = springs[i].getSpringStartID();  
        int springEndID = springs[i].getSpringEndID();  

        float springCoef  = springs[i].getSpringCoef();
        float damperCoef  = springs[i].getDamperCoef();
        float restLength  = springs[i].getSpringRestLength();

        Eigen::Vector3f posS = particles[springStartID].getPosition();
        Eigen::Vector3f posE = particles[springEndID].getPosition();
        Eigen::Vector3f velS = particles[springStartID].getVelocity();
        Eigen::Vector3f velE = particles[springEndID].getVelocity();

        Eigen::Vector3f sF = computeSpringForce(posS, posE, springCoef, restLength);
        Eigen::Vector3f dF = computeDamperForce(posS, posE, velS, velE, damperCoef);

        Eigen::Vector3f netF_S = sF + dF;
        Eigen::Vector3f netF_E = -netF_S;

        particles[springStartID].addForce(netF_S);
        particles[springEndID].addForce(netF_E);

    }
}


Eigen::Vector3f Jelly::computeSpringForce(const Eigen::Vector3f &positionA, const Eigen::Vector3f &positionB,
                                         const float springCoef, const float restLength) {
    // TODO#2-1: Compute spring force given the two positions of the spring.
    //   1. Review "particles.pptx" from p.9 - p.13
    //   2. The sample below just set spring force to zero
    Eigen::Vector3f Dir = positionB - positionA;
    Eigen::Vector3f sF = springCoef * (Dir.norm() - restLength) * (Dir / Dir.norm());
    return sF;
}

Eigen::Vector3f Jelly::computeDamperForce(const Eigen::Vector3f &positionA, const Eigen::Vector3f &positionB,
                                         const Eigen::Vector3f &velocityA, const Eigen::Vector3f &velocityB,
                                         const float damperCoef) {
    // TODO#2-2: Compute damper force given the two positions and the two velocities of the spring.
    //   1. Review "particles.pptx" from p.9 - p.13
    //   2. The sample below just set damper force to zero

    Eigen::Vector3f Dir = positionB - positionA;
    Eigen::Vector3f dir_N = Dir / Dir.norm();
    Eigen::Vector3f relativeVel = velocityB - velocityA;
    float speedAlongSpring = relativeVel.dot(dir_N);
    Eigen::Vector3f dF = damperCoef * speedAlongSpring * dir_N;

    return dF;
}


void Jelly::initializeParticle() {
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                Particle Particle;
                float offset_x = (float)((i - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));
                float offset_y = (float)((j - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));
                float offset_z = (float)((k - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));
                Particle.setPosition(Eigen::Vector3f(initialPosition(0) + offset_x, initialPosition(1) + offset_y,
                                                     initialPosition(2) + offset_z));
                particles.push_back(Particle);

            }
        }
    }
}

void Jelly::initializeSpring()
{
    // TODO#1: Connect particles with springs.  
    //   1. Consider the type of springs and compute indices of particles which the spring connect to.
    //   2. Compute rest spring length using particle positions.
    //   2. Iterate the particles. Push spring objects into `springs` vector
    // Note:
    //   1. The particles index can be computed in a similar way as below:
    //   ===============================================
    //   0 1 2 3 ... particlesPerEdge
    //   particlesPerEdge + 1 ....
    //   ... ... particlesPerEdge * particlesPerEdge - 1
    //   ===============================================
    // Here is a simple example which connects the structrual springs along z-axis.

    // struct
    // z-direction
    int struct_num = 0;
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge - 1; k++) {
              
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID = iParticleID + 1; 
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                Eigen::Vector3f SpringEndPos   = particles[iNeighborID].getPosition();
                Eigen::Vector3f Length         = SpringStartPos - SpringEndPos;
                float absLength = Length.norm();
                springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefStruct, damperCoefStruct,
                                         Spring::SpringType::STRUCT));
                struct_num++;
            }
        }
    }
  

    // x-direction
    for (int i = 0; i < particleNumPerEdge - 1; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID = (i + 1) * particleNumPerFace + j * particleNumPerEdge + k;
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                Eigen::Vector3f SpringEndPos   = particles[iNeighborID].getPosition();
                Eigen::Vector3f Length         = SpringStartPos - SpringEndPos;
                float absLength = Length.norm();

                springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefStruct, damperCoefStruct,
                                         Spring::SpringType::STRUCT));
                struct_num++;
            }
        }
    }

    // y-direction
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge - 1; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
              
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID = i * particleNumPerFace + (j + 1) * particleNumPerEdge + k;
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                Eigen::Vector3f SpringEndPos   = particles[iNeighborID].getPosition();
                Eigen::Vector3f Length         = SpringStartPos - SpringEndPos;
                float absLength = Length.norm();

                springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefStruct, damperCoefStruct,
                                         Spring::SpringType::STRUCT));
                struct_num++;
            }
        }
}                                                                                                                                                                                                                                           

    // bend
    int bend_num = 0;

    
    for (int i = 0; i < particleNumPerEdge - 2; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
              
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID = (i + 2) * particleNumPerFace + j * particleNumPerEdge + k;
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                Eigen::Vector3f SpringEndPos   = particles[iNeighborID].getPosition();
                Eigen::Vector3f Length         = SpringStartPos - SpringEndPos;
                float absLength = Length.norm();

                springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefBending, damperCoefBending,
                                         Spring::SpringType::BENDING));
                bend_num++;
            }
        }
    }

    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge - 2; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
              
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID = i * particleNumPerFace + (j + 2) * particleNumPerEdge + k;
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                Eigen::Vector3f SpringEndPos   = particles[iNeighborID].getPosition();
                Eigen::Vector3f Length         = SpringStartPos - SpringEndPos;
                float absLength = Length.norm();

                springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefBending, damperCoefBending,
                                         Spring::SpringType::BENDING));
                bend_num++;
            }
        }
    }

    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge - 2; k++) {
              
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID = iParticleID + 2; 
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                Eigen::Vector3f SpringEndPos   = particles[iNeighborID].getPosition();
                Eigen::Vector3f Length         = SpringStartPos - SpringEndPos;
                float absLength = Length.norm();

                springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefBending, damperCoefBending,
                                         Spring::SpringType::BENDING));
                bend_num++;
            }
        }
    }

    // shear
    int shear_num = 0;

    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                if (i < particleNumPerEdge - 1 && j < particleNumPerEdge - 1) {
                    int iNeighborID = (i + 1) * particleNumPerFace + (j + 1) * particleNumPerEdge + k;
                    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                    float absLength = (SpringStartPos - SpringEndPos).norm();
                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear, damperCoefShear,
                                             Spring::SpringType::SHEAR));
                    shear_num++;
                }
  
                if (i < particleNumPerEdge - 1 && j > 0) {
                    int iNeighborID = (i + 1) * particleNumPerFace + (j - 1) * particleNumPerEdge + k;
                    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                    float absLength = (SpringStartPos - SpringEndPos).norm();
                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear, damperCoefShear,
                                             Spring::SpringType::SHEAR));
                    shear_num++;
                }
            }
        }
    }

    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                if (i < particleNumPerEdge - 1 && k < particleNumPerEdge - 1) {
                    int iNeighborID = (i + 1) * particleNumPerFace + j * particleNumPerEdge + (k + 1);
                    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                    float absLength = (SpringStartPos - SpringEndPos).norm();
                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear, damperCoefShear,
                                             Spring::SpringType::SHEAR));
                    shear_num++;
                }
                if (i < particleNumPerEdge - 1 && k > 0) {
                    int iNeighborID = (i + 1) * particleNumPerFace + j * particleNumPerEdge + (k - 1);
                    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                    float absLength = (SpringStartPos - SpringEndPos).norm();
                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear, damperCoefShear,
                                             Spring::SpringType::SHEAR));
                    shear_num++;
                }
            }
        }
    }

    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                if (j < particleNumPerEdge - 1 && k < particleNumPerEdge - 1) {
                    int iNeighborID = i * particleNumPerFace + (j + 1) * particleNumPerEdge + (k + 1);
                    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                    float absLength = (SpringStartPos - SpringEndPos).norm();
                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear, damperCoefShear,
                                             Spring::SpringType::SHEAR));
                    shear_num++;
                }
                if (j < particleNumPerEdge - 1 && k > 0) {
                    int iNeighborID = i * particleNumPerFace + (j + 1) * particleNumPerEdge + (k - 1);
                    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                    float absLength = (SpringStartPos - SpringEndPos).norm();
                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear, damperCoefShear,
                                             Spring::SpringType::SHEAR));
                    shear_num++;
                }
            }
        }
    }

    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                if (i < particleNumPerEdge - 1 && j < particleNumPerEdge - 1 && k < particleNumPerEdge - 1) {
                    int iNeighborID = (i + 1) * particleNumPerFace
                                    + (j + 1) * particleNumPerEdge
                                    + (k + 1);
                    float absLength = (SpringStartPos - particles[iNeighborID].getPosition()).norm();
                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear, damperCoefShear,
                                             Spring::SpringType::SHEAR));
                    shear_num++;
                }
                if (i < particleNumPerEdge - 1 && j < particleNumPerEdge - 1 && k > 0) {
                    int iNeighborID = (i + 1) * particleNumPerFace
                                    + (j + 1) * particleNumPerEdge
                                    + (k - 1);
                    float absLength = (SpringStartPos - particles[iNeighborID].getPosition()).norm();
                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear, damperCoefShear,
                                             Spring::SpringType::SHEAR));
                    shear_num++;
                }

                if (i < particleNumPerEdge - 1 && j > 0 && k < particleNumPerEdge - 1) {
                    int iNeighborID = (i + 1) * particleNumPerFace
                                    + (j - 1) * particleNumPerEdge
                                    + (k + 1);
                    float absLength = (SpringStartPos - particles[iNeighborID].getPosition()).norm();
                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear, damperCoefShear,
                                             Spring::SpringType::SHEAR));
                    shear_num++;
                }
                if (i < particleNumPerEdge - 1 && j > 0 && k > 0) {
                    int iNeighborID = (i + 1) * particleNumPerFace
                                    + (j - 1) * particleNumPerEdge
                                    + (k - 1);
                    float absLength = (SpringStartPos - particles[iNeighborID].getPosition()).norm();
                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear, damperCoefShear,
                                             Spring::SpringType::SHEAR));
                    shear_num++;
                }
            }
        }
    }
}



void Jelly::updateSpringCoef(const float a_cdSpringCoef, const Spring::SpringType a_cSpringType) {
    for (unsigned int uiI = 0; uiI < springs.size(); uiI++) {
        if (springs[uiI].getType() == a_cSpringType) {
            springs[uiI].setSpringCoef(a_cdSpringCoef);
        }
    }
}

void Jelly::updateDamperCoef(const float a_cdDamperCoef, const Spring::SpringType a_cSpringType) {
    for (unsigned int uiI = 0; uiI < springs.size(); uiI++) {
        if (springs[uiI].getType() == a_cSpringType) {
            springs[uiI].setDamperCoef(a_cdDamperCoef);
        }
    }
}
}  //  namespace simulation
