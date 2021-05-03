#include "integrator.h"

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
    // TODO
    
    for (int i = 0; i < particleSystem.getCubeCount(); i++) {
        for (int j = 0; j < particleSystem.getCubePointer(i)->getParticleNum(); j++) {
            Particle currentParticle = particleSystem.getCubePointer(i)->getParticle(j);

            Eigen::Vector3f deltaPosition = currentParticle.getVelocity() * particleSystem.deltaTime;
            Eigen::Vector3f deltaVelocity = currentParticle.getAcceleration() * particleSystem.deltaTime;

            particleSystem.getCubePointer(i)->getParticle(j).addPosition(deltaPosition);
            particleSystem.getCubePointer(i)->getParticle(j).addVelocity(deltaVelocity);
        }
    }
}

//
// ImplicitEulerIntegrator
//

IntegratorType ImplicitEulerIntegrator::getType() { return IntegratorType::ImplicitEuler; }

void ImplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO

    for (int i = 0; i < particleSystem.getCubeCount(); i++) {
        Cube* currentCube = particleSystem.getCubePointer(i);

        particleSystem.getCubePointer(i);

        std::vector<Eigen::Vector3f> startingPosition;
        std::vector<Eigen::Vector3f> startingVelocity;
        Eigen::Vector3f deltaPosition;
        Eigen::Vector3f deltaVelocity;

        for (int j = 0; j < currentCube->getParticleNum(); j++) {
            Particle currentParticle = currentCube->getParticle(j);

            startingPosition.push_back(currentParticle.getPosition());
            startingVelocity.push_back(currentParticle.getVelocity());

            deltaPosition = currentParticle.getVelocity() * particleSystem.deltaTime;
            deltaVelocity = currentParticle.getAcceleration() * particleSystem.deltaTime;

            currentParticle.setPosition(startingPosition[j] + deltaPosition);
            currentParticle.setVelocity(startingVelocity[j] + deltaVelocity);
        }

        particleSystem.computeCubeForce(*currentCube);

        for (int j = 0; j < currentCube->getParticleNum(); j++) {
            Particle currentParticle = currentCube->getParticle(j);

            deltaPosition = currentParticle.getVelocity() * particleSystem.deltaTime;
            deltaVelocity = currentParticle.getAcceleration() * particleSystem.deltaTime;

            particleSystem.getCubePointer(i)->getParticle(j).setPosition(startingPosition[j] + deltaPosition);
            particleSystem.getCubePointer(i)->getParticle(j).setVelocity(startingVelocity[j] + deltaVelocity);
        }
    }
}

//
// MidpointEulerIntegrator
//

IntegratorType MidpointEulerIntegrator::getType() { return IntegratorType::MidpointEuler; }

void MidpointEulerIntegrator::integrate(MassSpringSystem& particleSystem) {

    // TODO
    // For midpoint euler, the deltaTime passed in is correct.
    // But this deltaTime is for a full step.
    // So you may need to adjust it before computing, but don't forget to restore original value.

    for (int i = 0; i < particleSystem.getCubeCount(); i++) {

        Cube* currentCube = particleSystem.getCubePointer(i);
        
        particleSystem.getCubePointer(i);

        std::vector<Eigen::Vector3f> startingPosition;
        std::vector<Eigen::Vector3f> startingVelocity;
        Eigen::Vector3f deltaPosition;
        Eigen::Vector3f deltaVelocity;

        for (int j = 0; j < currentCube->getParticleNum(); j++) {

            Particle currentParticle = currentCube->getParticle(j);

            startingPosition.push_back(currentParticle.getPosition());
            startingVelocity.push_back(currentParticle.getVelocity());

            deltaPosition = currentParticle.getVelocity() * particleSystem.deltaTime;
            deltaVelocity = currentParticle.getAcceleration() * particleSystem.deltaTime;

            currentParticle.setPosition(startingPosition[j] + deltaPosition / 2);
            currentParticle.setVelocity(startingVelocity[j] + deltaVelocity / 2);
        }

        particleSystem.computeCubeForce(*currentCube);

        for (int j = 0; j < currentCube->getParticleNum(); j++) {

            Particle currentParticle = currentCube->getParticle(j);

            deltaPosition = currentParticle.getVelocity() * particleSystem.deltaTime;
            deltaVelocity = currentParticle.getAcceleration() * particleSystem.deltaTime;

            particleSystem.getCubePointer(i)->getParticle(j).setPosition(startingPosition[j] + deltaPosition);
            particleSystem.getCubePointer(i)->getParticle(j).setVelocity(startingVelocity[j] + deltaVelocity);
        }
    }
}

//
// RungeKuttaFourthIntegrator
//

IntegratorType RungeKuttaFourthIntegrator::getType() { return IntegratorType::RungeKuttaFourth; }

void RungeKuttaFourthIntegrator::integrate(MassSpringSystem& particleSystem) {

    struct StateStep {
        Eigen::Vector3f deltaVel;
        Eigen::Vector3f deltaPos;
    };

    // TODO
    // StateStep struct is just a hint, you can use whatever you want.

    for (int i = 0; i < particleSystem.getCubeCount(); i++) {

        Cube* currentCube = particleSystem.getCubePointer(i);

        std::vector<Eigen::Vector3f> startingPosition;
        std::vector<Eigen::Vector3f> startingVelocity;
       
        std::vector<StateStep> k1;
        std::vector<StateStep> k2;
        std::vector<StateStep> k3;
        std::vector<StateStep> k4;

        for (int j = 0; j < currentCube->getParticleNum(); j++) {
            
            Particle currentParticle = currentCube->getParticle(j);
            StateStep tempStateStep;

            startingPosition.push_back(currentParticle.getPosition());
            startingVelocity.push_back(currentParticle.getVelocity());

            tempStateStep.deltaPos = currentParticle.getVelocity() * particleSystem.deltaTime;
            tempStateStep.deltaVel = currentParticle.getAcceleration() * particleSystem.deltaTime;

            k1.push_back(tempStateStep);

            currentParticle.setPosition(startingPosition[j] + k1[j].deltaPos / 2);
            currentParticle.setVelocity(startingVelocity[j] + k1[j].deltaVel / 2);
        }

        particleSystem.computeCubeForce(*currentCube);

        for (int j = 0; j < currentCube->getParticleNum(); j++) {

            Particle currentParticle = currentCube->getParticle(j);
            StateStep tempStateStep;

            tempStateStep.deltaPos = currentParticle.getVelocity() * particleSystem.deltaTime;
            tempStateStep.deltaVel = currentParticle.getAcceleration() * particleSystem.deltaTime;

            k2.push_back(tempStateStep);

            currentParticle.setPosition(startingPosition[j] + k2[j].deltaPos / 2);
            currentParticle.setVelocity(startingVelocity[j] + k2[j].deltaVel / 2);
        }
        particleSystem.computeCubeForce(*currentCube);

        for (int j = 0; j < currentCube->getParticleNum(); j++) {
            Particle currentParticle = currentCube->getParticle(j);
            StateStep tempStateStep;
            
            tempStateStep.deltaPos = currentParticle.getVelocity() * particleSystem.deltaTime;
            tempStateStep.deltaVel = currentParticle.getAcceleration() * particleSystem.deltaTime;

            k3.push_back(tempStateStep);

            currentParticle.setPosition(startingPosition[j] + k3[j].deltaPos);
            currentParticle.setVelocity(startingVelocity[j] + k3[j].deltaVel);
        }

        
        particleSystem.computeCubeForce(*currentCube);

        for (int j = 0; j < currentCube->getParticleNum(); j++) {
            Particle currentParticle = currentCube->getParticle(j);
            StateStep tempStateStep;

            tempStateStep.deltaPos = currentParticle.getVelocity() * particleSystem.deltaTime;
            tempStateStep.deltaVel = currentParticle.getAcceleration() * particleSystem.deltaTime;

            k4.push_back(tempStateStep);
        }

        for (int j = 0; j < currentCube->getParticleNum(); j++) {
            Particle currentParticle = currentCube->getParticle(j);

            Eigen::Vector3f deltaPosition =
                (k1[j].deltaPos + 2 * k2[j].deltaPos + 2 * k3[j].deltaPos + k4[j].deltaPos) / 6;
            Eigen::Vector3f deltaVelocity =
                (k1[j].deltaVel + 2 * k2[j].deltaVel + 2 * k3[j].deltaVel + k4[j].deltaVel) / 6;
            particleSystem.getCubePointer(i)->getParticle(j).setPosition(startingPosition[j] + deltaPosition);
            particleSystem.getCubePointer(i)->getParticle(j).setVelocity(startingVelocity[j] + deltaVelocity);
        }

        
    }

}
}  // namespace simulation
