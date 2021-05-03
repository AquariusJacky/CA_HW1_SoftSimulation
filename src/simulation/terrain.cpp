#include "terrain.h"

#include <stdexcept>
#include <cmath>

#include "../util/helper.h"

namespace simulation {
// Factory
std::unique_ptr<Terrain> TerrainFactory::CreateTerrain(TerrainType type) {
    switch (type) {
        case simulation::TerrainType::Plane:
            return std::make_unique<PlaneTerrain>();
        case simulation::TerrainType::Sphere:
            return std::make_unique<SphereTerrain>();
        case simulation::TerrainType::Bowl:
            return std::make_unique<BowlTerrain>();
        case simulation::TerrainType::TiltedPlane:
            return std::make_unique<TiltedPlaneTerrain>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}
// Terrain

Eigen::Matrix4f Terrain::getModelMatrix() { return modelMatrix; }

// Note:
// You should update each particles' velocity (base on the equation in
// slide) and force (contact force : resist + friction) in handleCollision function

// PlaneTerrain //
PlaneTerrain::PlaneTerrain() { modelMatrix = util::translate(0.0f, position[1], 0.0f) * util::scale(60, 1, 60); }

TerrainType PlaneTerrain::getType() { return TerrainType::Plane; }

void PlaneTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float resistCoef = 0.8f;
    constexpr float frictionCoef = 0.3f;
    
    // TODO
    
    for (int i = 0; i < cube.getParticleNum(); i++) {

        Particle currentParticle = cube.getParticle(i);
        
        if ((currentParticle.getPosition() - position).dot(normal) < eEPSILON) {

            float normalLength = normal.norm();

            Eigen::Vector3f totalVelocity = currentParticle.getVelocity();
            Eigen::Vector3f normalVelocity = (totalVelocity.dot(normal) / (normalLength * normalLength)) * normal;
            Eigen::Vector3f tangentialVelocity = totalVelocity - normalVelocity;

            if (currentParticle.getVelocity().dot(normal) < 0) {
                Eigen::Vector3f newVelocity = -resistCoef * normalVelocity + tangentialVelocity;
                cube.getParticle(i).setVelocity(newVelocity);
            }

            Eigen::Vector3f currentForce = currentParticle.getForce();

            if (currentForce.dot(normal) < 0) {
                Eigen::Vector3f collisionForce = -currentForce.dot(normal) * normal;
                Eigen::Vector3f frictionalForce =
                    frictionCoef * currentForce.dot(normal) * (tangentialVelocity / tangentialVelocity.norm());

                cube.getParticle(i).addForce(collisionForce);
                cube.getParticle(i).addForce(frictionalForce);
            }
        }
    }
}

// SphereTerrain //
SphereTerrain::SphereTerrain() { modelMatrix = util::translate(position) * util::scale(radius, radius, radius); }

TerrainType SphereTerrain::getType() { return TerrainType::Sphere; }

void SphereTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float resistCoef = 0.8f;
    constexpr float frictionCoef = 0.3f;

    // TODO

    for (int i = 0; i < cube.getParticleNum(); i++) {
        Particle currentParticle = cube.getParticle(i);
        Eigen::Vector3f terrainNormalVector = currentParticle.getPosition() - position;
        float centerToParticleLength = terrainNormalVector.norm();
        terrainNormalVector = terrainNormalVector / centerToParticleLength;

        if ((centerToParticleLength - radius) < eEPSILON) {

            Eigen::Vector3f totalVelocity = currentParticle.getVelocity();
            Eigen::Vector3f normalVelocity =
                (totalVelocity.dot(terrainNormalVector) / (centerToParticleLength * centerToParticleLength)) *
                (terrainNormalVector);
            Eigen::Vector3f tangentialVelocity = totalVelocity - normalVelocity;

            if (currentParticle.getVelocity().dot(terrainNormalVector) < 0) {
                float particleMass = currentParticle.getMass();

                Eigen::Vector3f newVelocity =
                    ((particleMass - mass) / (particleMass + mass)) * normalVelocity + tangentialVelocity;
                cube.getParticle(i).setVelocity(newVelocity);
            }

            Eigen::Vector3f currentForce = currentParticle.getForce();

            if (currentForce.dot(terrainNormalVector) < 0) {
                Eigen::Vector3f collisionForce = -currentForce.dot(terrainNormalVector) * terrainNormalVector;
                Eigen::Vector3f frictionalForce = frictionCoef * currentForce.dot(terrainNormalVector) *
                                                  (tangentialVelocity / tangentialVelocity.norm());

                cube.getParticle(i).addForce(collisionForce);
                cube.getParticle(i).addForce(frictionalForce);
            }
        }
    }

}

// BowlTerrain //

BowlTerrain::BowlTerrain() { modelMatrix = util::translate(position) * util::scale(radius, radius, radius); }

TerrainType BowlTerrain::getType() { return TerrainType::Bowl; }

void BowlTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float resistCoef = 0.8f;
    constexpr float frictionCoef = 0.3f;
    // TODO

    for (int i = 0; i < cube.getParticleNum(); i++) {
        Particle currentParticle = cube.getParticle(i);
        Eigen::Vector3f terrainNormalVector = position - currentParticle.getPosition();
        float centerToParticleLength = terrainNormalVector.norm();
        terrainNormalVector = terrainNormalVector / centerToParticleLength;

        if ((radius - centerToParticleLength) < eEPSILON) {

            Eigen::Vector3f totalVelocity = currentParticle.getVelocity();
            Eigen::Vector3f normalVelocity =
                (totalVelocity.dot(terrainNormalVector) / (centerToParticleLength * centerToParticleLength)) *
                (terrainNormalVector);
            Eigen::Vector3f tangentialVelocity = totalVelocity - normalVelocity;

            if (totalVelocity.dot(terrainNormalVector) < 0) {

                float particleMass = currentParticle.getMass();

                Eigen::Vector3f newVelocity = 
                    ((particleMass - mass) / (particleMass + mass)) * normalVelocity + tangentialVelocity;
                cube.getParticle(i).setVelocity(newVelocity);
            }

            Eigen::Vector3f currentForce = currentParticle.getForce();
            
            if (currentForce.dot(terrainNormalVector) < 0) {
                Eigen::Vector3f collisionForce = -currentForce.dot(terrainNormalVector) * terrainNormalVector;
                Eigen::Vector3f frictionalForce = frictionCoef * currentForce.dot(terrainNormalVector) *
                                                  (tangentialVelocity / tangentialVelocity.norm());

                cube.getParticle(i).addForce(collisionForce);
                cube.getParticle(i).addForce(frictionalForce);
            }
        }
    }
}

// TiltedPlaneTerrain //

TiltedPlaneTerrain::TiltedPlaneTerrain() { modelMatrix = util::rotateDegree(0, 0, -45) * util::scale(60, 1, 60); }

TerrainType TiltedPlaneTerrain::getType() { return TerrainType::TiltedPlane; }

void TiltedPlaneTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float resistCoef = 0.8f;
    constexpr float frictionCoef = 0.3f;

    // TODO

    for (int i = 0; i < cube.getParticleNum(); i++) {
        Particle currentParticle = cube.getParticle(i);

        if ((currentParticle.getPosition() - position).dot(normal) < eEPSILON) {
            float normalLength = normal.norm();

            Eigen::Vector3f totalVelocity = currentParticle.getVelocity();
            Eigen::Vector3f normalVelocity = (totalVelocity.dot(normal) / (normalLength * normalLength)) * normal;
            Eigen::Vector3f tangentialVelocity = totalVelocity - normalVelocity;

            if (currentParticle.getVelocity().dot(normal) < 0) {
                Eigen::Vector3f newVelocity = -resistCoef * normalVelocity + tangentialVelocity;
                cube.getParticle(i).setVelocity(newVelocity);
            }

            Eigen::Vector3f currentForce = currentParticle.getForce();

            if (currentForce.dot(normal) < 0) {
                Eigen::Vector3f collisionForce = -currentForce.dot(normal) * normal;
                Eigen::Vector3f frictionalForce =
                    frictionCoef * currentForce.dot(normal) * (tangentialVelocity / tangentialVelocity.norm());

                cube.getParticle(i).addForce(collisionForce);
                cube.getParticle(i).addForce(frictionalForce);
            }
        }
    }
}
}  // namespace simulation
