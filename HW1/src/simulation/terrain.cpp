#include "terrain.h"
#include <stdexcept>
#include "../util/helper.h"
#include <iostream>
#include <cmath>  

#define M_PI 3.14159265358979323846
namespace simulation {
// Factory
std::unique_ptr<Terrain> TerrainFactory::CreateTerrain(TerrainType type) {
    switch (type) {
        case simulation::TerrainType::Plane:
            return std::make_unique<PlaneTerrain>();
        case simulation::TerrainType::Elevator:
            return std::make_unique<ElevatorTerrain>();

        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}
// Terrain

Eigen::Matrix4f Terrain::getModelMatrix() { return modelMatrix; }

float Terrain::getMass() const { return mass; }

Eigen::Vector3f Terrain::getPosition() const { return position; }

Eigen::Vector3f Terrain::getVelocity() const { return velocity; }

Eigen::Vector3f Terrain::getAcceleration() const { return force / mass; }

Eigen::Vector3f Terrain::getForce() const { return force; }

void Terrain::setMass(const float _mass) { mass = _mass; }

void Terrain::setPosition(const Eigen::Vector3f& _position) { 
    modelMatrix = util::translate(_position - position) * modelMatrix;
    position = _position;
}

void Terrain::setVelocity(const Eigen::Vector3f& _velocity) { velocity = _velocity; }

void Terrain::setAcceleration(const Eigen::Vector3f& _acceleration) { force = _acceleration * mass; }

void Terrain::setForce(const Eigen::Vector3f& _force) { force = _force; }

void Terrain::addPosition(const Eigen::Vector3f& _position) { 
    position += _position;
    modelMatrix = util::translate(_position) * modelMatrix;
}

void Terrain::addVelocity(const Eigen::Vector3f& _velocity) { velocity += _velocity; }

void Terrain::addAcceleration(const Eigen::Vector3f& _acceleration) { force += _acceleration * mass; }

void Terrain::addForce(const Eigen::Vector3f& _force) { force += _force; }

// Note:
// You should update each particles' velocity (base on the equation in
// slide) and force (contact force : resist + friction) in handleCollision function

// PlaneTerrain //

PlaneTerrain::PlaneTerrain() { reset(); }

void PlaneTerrain::reset() { 
    modelMatrix = util::translate(0.0f, position[1], 0.0f) * util::rotateDegree(0, 0, -20) * util::scale(30, 1, 30);
}


TerrainType PlaneTerrain::getType() { return TerrainType::Plane; }

void PlaneTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.5f;
    // TODO#3-1: Handle collision when a particle collide with the plane terrain.
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data.
    //   3. The plane spans 30x30 units in the XZ plane and is rotated -20 degrees around the Z-axis.
    // Hint:
    //   1. Review "particles.pptx" from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.
    const float angle = 20.0f * static_cast<float>(M_PI) / 180.0f; 
    Eigen::Vector3f planeN(std::sin(angle), std::cos(angle), 0.0f);
    
    for (int i = 0; i < jelly.getParticleNum(); i++) {
        Particle& particle = jelly.getParticle(i);
        Eigen::Vector3f pos = particle.getPosition();
        
        if (abs(pos.x()) > 15.0f || abs(pos.z()) > 15.0f ) continue; 
        float dist = pos.dot(planeN);
        if (dist < eEPSILON) {
            float penetration = eEPSILON - dist;

            Eigen::Vector3f Vel = particle.getVelocity();
            float vDot = Vel.dot(planeN);
            Eigen::Vector3f VelN = vDot * planeN;          
            Eigen::Vector3f VelT = Vel - VelN;     

            if (vDot < 0.0f) {
                VelN = -coefResist * VelN;
            }
            VelT = (1.0f - coefFriction) * VelT;
            particle.setVelocity(VelN + VelT);
            Eigen::Vector3f contactF = (penetration / delta_T) * planeN;
            particle.addForce(contactF);
        }
    }
}

ElevatorTerrain::ElevatorTerrain() { 
    reset();
}

void ElevatorTerrain::reset() {
    modelMatrix = util::translate(0.0f, 1.0f, 0.0f) * util::rotateDegree(0, 0, 0) * util::scale(5, 1, 5);
    position = Eigen::Vector3f(0.0f, 1.0f, 0.0f);
    velocity = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
}

TerrainType ElevatorTerrain::getType() { return TerrainType::Elevator; }

void ElevatorTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.5f;

    // TODO#3-2: Implement the collision handling between the jelly and the elevator
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data.
    //   3. The elevator plane spans 5x5 units in the XZ plane.
    // Hint:
    //   1. Review "particles.pptx" from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.
    Eigen::Vector3f planeN(0.0f, 1.0f, 0.0f);
    Eigen::Vector3f elevatorPos = this->getPosition();
    for (int i = 0; i < jelly.getParticleNum(); i++) {
        Particle& particle = jelly.getParticle(i);
        Eigen::Vector3f pos = particle.getPosition();
        if (abs(pos.x() - elevatorPos.x())  < -2.5f || abs(pos.z() - elevatorPos.z()) < -2.5f ) continue;  

        float dist = planeN.dot(pos - elevatorPos);
        Eigen::Vector3f elevatorVel = this->getVelocity();
        Eigen::Vector3f Vel = particle.getVelocity() - elevatorVel;
        float vDot = planeN.dot(Vel);  
        if (dist < -eEPSILON && vDot < 0.0f) {
            Eigen::Vector3f VelN = vDot * planeN;          
            Eigen::Vector3f VelT = Vel - VelN;
            particle.setVelocity(VelT - (coefResist)*VelN + elevatorVel);
            
            float fvDot = planeN.dot(particle.getForce());

            if (fvDot < 0.0f) {
                Eigen::Vector3f contactF = -fvDot * planeN;
                particle.addForce(contactF);
                Eigen::Vector3f tVel = (Vel- elevatorVel) - planeN * vDot;
                Eigen::Vector3f frictionDir = tVel.normalized();
                Eigen::Vector3f frictionF = -coefFriction * (-fvDot) * frictionDir;
                particle.addForce(frictionF);

            }
        }
    }
}

}  // namespace simulation
