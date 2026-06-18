#include "integrator.h"
#include <iostream>

#include <vector>
namespace simulation {
// Factory
std::unique_ptr<Integrator> IntegratorFactory::CreateIntegrator(IntegratorType type) {
    switch (type) {
        case simulation::IntegratorType::ExplicitEuler:
            return std::make_unique<ExplicitEulerIntegrator>();
        case simulation::IntegratorType::ImplicitEuler:
            return std::make_unique<ImplicitEulerIntegrator>();
        case simulation::IntegratorType::MidpointEuler:
            return std::make_unique<MidpointEulerIntegrator>();
        case simulation::IntegratorType::RungeKuttaFourth:
            return std::make_unique<RungeKuttaFourthIntegrator>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}

//
// ExplicitEulerIntegrator
//

IntegratorType ExplicitEulerIntegrator::getType() { return IntegratorType::ExplicitEuler; }

void ExplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-1: Integrate position and velocity
    //   1. Integrate position using current velocity.
    //   2. Integrate velocity using current acceleration.
    //   3. Clear force
    // Note:
    //   1. You should do this first because it is very simple. Then you can check whether your collision is correct or
    //   not.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Review "ODE_basics.pptx" from p.15 - p.16
    double dt = particleSystem.deltaTime;
    for (auto &jelly : particleSystem.jellies) {
        int n = jelly.getParticleNum();
        for (int i = 0; i < n; ++i) {
            Particle &p = jelly.getParticle(i);
            p.setPosition(p.getPosition() + dt * p.getVelocity());
            p.setVelocity(p.getVelocity() + dt * p.getAcceleration());
            p.setForce(Eigen::Vector3f::Zero());
        }
    }

    Terrain* elevator = particleSystem.elevatorTerrain.get();
    if (elevator) {
        elevator->setPosition(elevator->getPosition() + dt * elevator->getVelocity());
        elevator->setVelocity(elevator->getVelocity() + dt * elevator->getAcceleration());
        elevator->setForce(Eigen::Vector3f::Zero());
    }
}

//
// ImplicitEulerIntegrator
//

IntegratorType ImplicitEulerIntegrator::getType() { return IntegratorType::ImplicitEuler; }

void ImplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-2: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 and Vn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce`and `MassSpringSystem::computeElevatorForce`
    //      with modified position and velocity to get Xn+1.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Remember to reset particleSystem.elevatorCounter back to the original state
    //   5. Review "ODE_implicit.pptx" from p.18 - p.19
    double dt = particleSystem.deltaTime;
    int originalElevatorCounter = particleSystem.elevatorCounter;
    std::vector<Eigen::Vector3f> backupPos, backupVel;
    Eigen::Vector3f elevatorBackupPos = particleSystem.elevatorTerrain->getPosition();
    Eigen::Vector3f elevatorBackupVel = particleSystem.elevatorTerrain->getVelocity();
    for (auto &jelly : particleSystem.jellies) {
        int n = jelly.getParticleNum();
        for (int i = 0; i < n; ++i) {
            Particle &p = jelly.getParticle(i);
            backupPos.push_back(p.getPosition());
            backupVel.push_back(p.getVelocity());
        }

        for (int i = 0; i < n; ++i) {
            Particle &p = jelly.getParticle(i);
            p.setPosition(p.getPosition() + dt * p.getVelocity());
            p.setVelocity(p.getVelocity() + dt * p.getAcceleration());
            p.setForce(Eigen::Vector3f::Zero());
        }
    }

    Terrain* elevator = particleSystem.elevatorTerrain.get();
    if (elevator) {
        elevator->setVelocity(elevator->getVelocity() + dt * elevator->getAcceleration());
        elevator->setPosition(elevator->getPosition() + dt * elevator->getVelocity());
        elevator->setForce(Eigen::Vector3f::Zero());
    }
    particleSystem.computeAllForce();
    particleSystem.computeElevatorForce();
    int j = 0 ;
    for (auto &jelly : particleSystem.jellies) {
        int n = jelly.getParticleNum();
        for (int i = 0; i < n; ++i) {
            Particle &p = jelly.getParticle(i);
            p.setVelocity(backupVel[j] + dt * p.getAcceleration());
            p.setPosition(backupPos[j++] + dt * p.getVelocity());
            p.setForce(Eigen::Vector3f::Zero());
        }
    }
    
    if (elevator) {
        elevator->setPosition(elevatorBackupPos + dt * elevator->getVelocity());
        elevator->setVelocity(elevatorBackupVel + dt * elevator->getAcceleration());
        elevator->setForce(Eigen::Vector3f::Zero());
    }
    particleSystem.elevatorCounter = originalElevatorCounter;
}

//
// MidpointEulerIntegrator
//

IntegratorType MidpointEulerIntegrator::getType() { return IntegratorType::MidpointEuler; }

void MidpointEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-3: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce`and `MassSpringSystem::computeElevatorForce`
    //      with modified position and velocity to get Xn+1.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Remember to reset particleSystem.elevatorCounter back to the original state
    //   5. Review "ODE_basics.pptx" from p.18 - p.19
    double dt = particleSystem.deltaTime ;
    int originalElevatorCounter = particleSystem.elevatorCounter;
    std::vector<Eigen::Vector3f> backupPos, backupVel, MidpointVel;
    for (auto &jelly : particleSystem.jellies) {
        int n = jelly.getParticleNum();
        for (int i = 0; i < n; ++i) {
            Particle &p = jelly.getParticle(i);
            backupPos.push_back(p.getPosition());
            backupVel.push_back(p.getVelocity());
        }
        for (int i = 0; i < n; ++i) {
            Particle &p = jelly.getParticle(i);
            p.setPosition(p.getPosition() + dt * 0.5f * p.getVelocity());
            MidpointVel.push_back(p.getVelocity() + dt * 0.5f * p.getAcceleration());
            p.setVelocity(p.getVelocity() + dt * 0.5f * p.getAcceleration());
            p.setForce(Eigen::Vector3f::Zero());
        }
    }
    
    Terrain* elevator = particleSystem.elevatorTerrain.get();
    Eigen::Vector3f elevatorBackupPos = particleSystem.elevatorTerrain->getPosition();
    Eigen::Vector3f elevatorBackupVel = particleSystem.elevatorTerrain->getVelocity();
    if (elevator) {
        elevator->setVelocity(elevatorBackupVel + dt * 0.5f * elevator->getAcceleration() );
        elevator->setPosition(elevatorBackupPos + dt * 0.5f * elevator->getVelocity());
        elevator->setForce(Eigen::Vector3f::Zero());
    }    
    particleSystem.computeAllForce();
    particleSystem.computeElevatorForce();
    int j = 0 ;
    
    for (auto &jelly : particleSystem.jellies) {
        int n = jelly.getParticleNum();
        for (int i = 0; i < n; ++i) {
            Particle &p = jelly.getParticle(i);
            p.setVelocity(backupVel[j] + dt * 0.5f * p.getAcceleration());
            p.setPosition(backupPos[j] + dt * 0.5f * MidpointVel[j++]);
            p.setForce(Eigen::Vector3f::Zero());
        }
    }


    if (elevator) {
        elevator->setPosition(elevatorBackupPos + dt * (elevatorBackupVel + dt * 0.5f * elevator->getAcceleration() ));
        elevator->setVelocity(elevatorBackupVel + dt * elevator->getAcceleration());
        elevator->setForce(Eigen::Vector3f::Zero());
    }
    particleSystem.elevatorCounter = originalElevatorCounter;
}

//
// RungeKuttaFourthIntegrator
//

IntegratorType RungeKuttaFourthIntegrator::getType() { return IntegratorType::RungeKuttaFourth; }

void RungeKuttaFourthIntegrator::integrate(MassSpringSystem& particleSystem) {
    struct StateStep {
        Eigen::Vector3f deltaVel = Eigen::Vector3f::Zero();
        Eigen::Vector3f deltaPos = Eigen::Vector3f::Zero();
    };
    // TODO#4-4: Integrate velocity and acceleration
    //   1. Backup original particles' data.
    //   2. Compute k1, k2, k3, k4
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce`and `MassSpringSystem::computeElevatorForce`
    //      with modified position and velocity to get Xn+1.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Remember to reset particleSystem.elevatorCounter back to the original state
    //   5. StateStep struct is just a hint, you can use whatever you want.
    //   6. Review "ODE_basics.pptx" from p.21
    double dt = particleSystem.deltaTime;
    int originalElevatorCounter = particleSystem.elevatorCounter;
    std::vector<Eigen::Vector3f> backupPos, backupVel;

    for (auto &jelly : particleSystem.jellies) {
        int n = jelly.getParticleNum();
        for (int i = 0; i < n; ++i) {
            Particle &p = jelly.getParticle(i);
            backupPos.push_back(p.getPosition());
            backupVel.push_back(p.getVelocity());
        }
    }
    Terrain* elevator = particleSystem.elevatorTerrain.get();
    Eigen::Vector3f elevatorBackupPos = particleSystem.elevatorTerrain->getPosition();
    Eigen::Vector3f elevatorBackupVel = particleSystem.elevatorTerrain->getVelocity();
    std::vector<StateStep> k1, k2, k3, k4;
    StateStep k1_e, k2_e, k3_e, k4_e ; 
    int j = 0 ; 
    for (auto &jelly : particleSystem.jellies) {
        int n = jelly.getParticleNum();
        for (int i = 0; i < n; ++i) {
            Particle &p = jelly.getParticle(i);
            p.setPosition(backupPos[j] + dt * 0.5f * p.getVelocity());
            p.setVelocity(backupVel[j++] + dt * 0.5f *p.getAcceleration());
            k1.push_back(StateStep{dt * p.getAcceleration(), dt * p.getVelocity()});
        }
    }
    if (elevator) {
        elevator->setPosition(elevatorBackupPos + dt * 0.5f * elevator->getVelocity());
        elevator->setVelocity(elevatorBackupVel + dt * 0.5f * elevator->getAcceleration());
        k1_e = StateStep{dt * elevator->getAcceleration(), dt * elevator->getVelocity()};
        elevator->setForce(Eigen::Vector3f::Zero());
    }
    particleSystem.computeAllForce();
    particleSystem.computeElevatorForce();
    
    j = 0 ;
    for (auto &jelly : particleSystem.jellies) {
        int n = jelly.getParticleNum();
        for (int i = 0; i < n; ++i) {
            Particle &p = jelly.getParticle(i);
            p.setPosition(backupPos[j] + dt * 0.5f * p.getVelocity());
            p.setVelocity(backupVel[j++] + dt * 0.5f *p.getAcceleration());
            k2.push_back(StateStep{dt * p.getAcceleration(), dt * p.getVelocity()});
        }
    }
    if (elevator) {
        elevator->setPosition(elevatorBackupPos + dt * 0.5f * elevator->getVelocity());
        elevator->setVelocity(elevatorBackupVel + dt * 0.5f * elevator->getAcceleration());
        k2_e = StateStep{dt * elevator->getAcceleration(), dt * elevator->getVelocity()};
        elevator->setForce(Eigen::Vector3f::Zero());
    }
    particleSystem.computeAllForce();
    particleSystem.computeElevatorForce();

    j = 0 ;
    for (auto &jelly : particleSystem.jellies) {
        int n = jelly.getParticleNum();
        for (int i = 0; i < n; ++i) {
            Particle &p = jelly.getParticle(i);
            p.setPosition(backupPos[j] + dt * 0.5f * p.getVelocity());
            p.setVelocity(backupVel[j++] + dt * 0.5f *p.getAcceleration());
            k3.push_back(StateStep{dt * p.getAcceleration(), dt * p.getVelocity()});
        }
    }
    if (elevator) {
        elevator->setPosition(elevatorBackupPos + dt * 0.5f * elevator->getVelocity());
        elevator->setVelocity(elevatorBackupVel + dt * 0.5f * elevator->getAcceleration());
        k3_e = StateStep{dt * elevator->getAcceleration(), dt * elevator->getVelocity()};
        elevator->setForce(Eigen::Vector3f::Zero());
    }
    particleSystem.computeAllForce();
    particleSystem.computeElevatorForce();

    for (auto &jelly : particleSystem.jellies) {
        int n = jelly.getParticleNum();
        for (int i = 0; i < n; ++i) {
            Particle &p = jelly.getParticle(i);
            k4.push_back(StateStep{dt * p.getAcceleration(), dt * p.getVelocity()});
        }
    }
    if (elevator) {
        k4_e = StateStep{dt * elevator->getAcceleration(), dt * elevator->getVelocity()};
    }

    j = 0 ; 
    for (auto &jelly : particleSystem.jellies) {
        int n = jelly.getParticleNum();
        for (int i = 0; i < n; ++i) {
            Particle &p = jelly.getParticle(i);
            p.setPosition(backupPos[j] + (k1[j].deltaPos + 2.0f * k2[j].deltaPos + 2.0f * k3[j].deltaPos + k4[j].deltaPos) / 6.0f);
            p.setVelocity(backupVel[j] + (k1[j].deltaVel + 2.0f * k2[j].deltaVel + 2.0f * k3[j].deltaVel + k4[j++].deltaVel) / 6.0f);
            p.setForce(Eigen::Vector3f::Zero());   
        }
    }
    if (elevator) {
        elevator->setPosition(elevatorBackupPos + (k1_e.deltaPos + 2.0f * k2_e.deltaPos + 2.0f * k3_e.deltaPos + k4_e.deltaPos) / 6.0f);
        elevator->setVelocity(elevatorBackupVel + (k1_e.deltaVel + 2.0f * k2_e.deltaVel + 2.0f * k3_e.deltaVel + k4_e.deltaVel) / 6.0);
        elevator->setForce(Eigen::Vector3f::Zero());
    }
    particleSystem.elevatorCounter = originalElevatorCounter;
}


}  // namespace simulation
