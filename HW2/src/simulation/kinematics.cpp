#include "simulation/kinematics.h"

#include <cmath>
#include <iostream>
#include "Eigen/Dense"
#include "acclaim/bone.h"
#include "util/helper.h"

namespace kinematics {

void forwardSolver(const acclaim::Posture& posture, acclaim::Bone* bone) {
    // TODO#1: Forward Kinematic
    // Hint:
    // - Traverse the skeleton tree from root to leaves.
    // - Compute each bone's global rotation and global position.
    // - Use local rotation (from posture) and bone hierarchy (parent rotation, offset, etc).
    // - Remember to update both bone->start_position and bone->end_position.
    // - Use bone->rotation to store global rotation (after combining parent, local, etc).
    Eigen::Affine3d localRot(util::rotateDegreeZYX(posture.bone_rotations[bone->idx]));

    if (bone->parent) {
        bone->start_position = bone->parent->end_position;
        bone->rotation = bone->parent->rotation * bone->rot_parent_current * localRot;
    } else {
        bone->start_position = posture.bone_translations[bone->idx];
        bone->rotation = bone->rot_parent_current * localRot;
    }

    Eigen::Vector3d dir = bone->dir.head<3>();
    Eigen::Vector3d offset = bone->rotation * (dir * bone->length);
    bone->end_position.head<3>() = bone->start_position.head<3>() + offset;
    bone->end_position[3] = 1.0;

    if (bone->child) forwardSolver(posture, bone->child);
    if (bone->sibling) forwardSolver(posture, bone->sibling);
}


Eigen::VectorXd pseudoInverseLinearSolver(const Eigen::Matrix4Xd& Jacobian, const Eigen::Vector4d& target) {
    Eigen::VectorXd deltatheta;
    // TODO#2: Inverse linear solver (find x which min(| jacobian * x - target |))
    // Hint:
    //   1. Linear algebra - least squares solution
    //   2. https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Construction
    // Note:
    //   1. SVD or other pseudo-inverse method is useful
    //   2. Some of them have some limitation, if you use that method you should check it.
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Jacobian, Eigen::ComputeThinU | Eigen::ComputeThinV);
    deltatheta = svd.solve(target);
    return deltatheta;
}

/**
 * @brief Perform inverse kinematics (IK)
 *
 * @param target_pos The position where `end_bone` will move to.
 * @param obs_pos The position where the obstacle is at
 * @param obsActive Whether the obstacle is active or not
 * @param start_bone This bone is the last bone you can move while doing IK
 * @param end_bone This bone will try to reach `target_pos`
 * @param posture The original AMC motion's reference, you need to modify this
 *
 * @return True if IK is stable (HW2 bonus)
 */
bool inverseJacobianIKSolver(const Eigen::Vector4d& target_pos, const Eigen::Vector4d& obs_pos, bool obsActive,
                             acclaim::Bone* start_bone, acclaim::Bone* end_bone, acclaim::Posture& posture) {
    constexpr int max_iteration = 1000;
    constexpr double epsilon = 1E-3;
    constexpr double step = 0.1;
    constexpr double obsAvoidThreshold = 1.01;  // if bone is within 1 unit from obstacle

    // Since bone stores in bones[i] that i == bone->idx, we can use bone - bone->idx to find bones[0] which is the
    // root.
    acclaim::Bone* root_bone = start_bone - start_bone->idx;

    // TODO#3:
    // Perform inverse kinematics (IK)
    // HINTs will tell you what should do in that area.
    // Of course you can ignore it (Any code below this line) and write your own code.
    acclaim::Posture original_posture(posture);

    // TODO#3-1:
    // Calculate number of bones need to move to perform IK, store in `bone_num`
    // (a.k.a. how may bones from end_bone to its parent than to start_bone (include both side))
    // Store the bones need to move to perform IK into boneList
    // Hint:
    //   1. Traverse from end_bone to start_bone is easier than start to end (since there is only 1 parent)
    //   2. If start bone is not reachable from end. Go to root first.
    // Note:
    //   1. Both start and end should be in the list

    size_t bone_num = 0;
    std::vector<acclaim::Bone *> boneList;
    for (acclaim::Bone *b = end_bone; b; b = b->parent) {
        boneList.push_back(b);
        if (b == start_bone) break;
    }
    bone_num = boneList.size();
    if (bone_num == 0) {
        return false;
    }
    acclaim::Bone *current = end_bone;
    Eigen::Matrix4Xd Jacobian(4, 3 * bone_num);
    Jacobian.setZero();
    bool limitExceeded = false; //bonus
    for (int iter = 0; iter < max_iteration; ++iter) {
        forwardSolver(posture, root_bone);
        Eigen::Vector4d desiredVector = target_pos - end_bone->end_position;

        if (desiredVector.head<3>().norm() < epsilon) {
            return true;
        }
        // TODO#3-2 (compute jacobian)
        //   1. Compute arm vectors
        //   2. Compute jacobian columns, store in `Jacobian`
        // Hint:
        //   1. You should not put rotation in jacobian if it doesn't have that DoF.
        //   2. jacobian.col(/* some column index */) = /* jacobian column */

        for (size_t i = 0; i < bone_num; ++i) {
            acclaim::Bone* b = boneList[i];
            Eigen::Vector3d jointPos = b->start_position.head<3>();
            Eigen::Matrix3d R = b->rotation.linear();
            
            if (b->dofrx) {
                Eigen::Vector3d u = R * Eigen::Vector3d::UnitX();
                Jacobian.block<3,1>(0, 3*i + 0) = u.cross(end_bone->end_position.head<3>() - jointPos);
            }
            if (b->dofry) {
                Eigen::Vector3d u = R * Eigen::Vector3d::UnitY();
                Jacobian.block<3,1>(0, 3*i + 1) = u.cross(end_bone->end_position.head<3>() - jointPos);
            }
            if (b->dofrz) {
                Eigen::Vector3d u = R * Eigen::Vector3d::UnitZ();
                Jacobian.block<3,1>(0, 3*i + 2) = u.cross(end_bone->end_position.head<3>() - jointPos);
            }
        }
        // TODO#3-3 (obstacle avoidance)
        //  1. Iterate through all bones in `boneList`.
        //  2. Compute the center of each bone (average of start and end positions).
        //  3. Calculate the vector from obstacle center to bone center.
        //  4. If distance is below threshold, compute repulsive vector.
        //  5. Add this repulsive vector to `desiredVector`.
        // Hint:
        // - Use a constant threshold distance to determine proximity.
        // - The repulsive vector should point away from the obstacle.
        // - Use `.head<3>().norm()` to compute 3D distance from a 4D vector.
        // - Normalize the repulsive vector and scale it based on how close it is.
        if (obsActive) {
            for (acclaim::Bone* b : boneList) {
                Eigen::Vector4d center = 0.5 * (b->start_position + b->end_position);
                Eigen::Vector3d obsDiff = center.head<3>() - obs_pos.head<3>();
                double dist = obsDiff.norm();
                if (dist < obsAvoidThreshold) {
                    if (desiredVector.head<3>().dot(obsDiff) < 0) {
                        Eigen::Vector3d repulse = obsDiff.normalized() * (obsAvoidThreshold - dist);
                        desiredVector.head<3>() += repulse;
                    }
                }
            }
        }

        Eigen::VectorXd deltatheta(3 * bone_num);
        deltatheta = step * pseudoInverseLinearSolver(Jacobian, desiredVector);
        // TODO#3-4 (update rotation)
        //   Update `posture.bone_rotation` (in euler angle / degrees) using deltaTheta
        // Hint:
        //   1. You can ignore rotation limit of the bone.
        // Bonus:
        //   1. You cannot ignore rotation limit of the bone.
        for (int i = 0; i < bone_num; ++i) {
            acclaim::Bone* b = boneList[i];
            auto &euler = posture.bone_rotations[b->idx];

            if (b->dofrx) {
                double before = euler[0];
                euler[0] += deltatheta(3 * i + 0) * 180.0 / std::acos(-1) ;
                euler[0] = std::clamp(euler[0], double(b->rxmin), double(b->rxmax));
                if (euler[0] != before) limitExceeded = true; //bonus
            }
            if (b->dofry) {
                double before = euler[1];
                euler[1] += deltatheta(3 * i + 1) * 180.0 / std::acos(-1) ;
                euler[1] = std::clamp(euler[1], double(b->rymin), double(b->rymax));
                if (euler[1] != before) limitExceeded = true; //bonus
            }
            if (b->dofrz) {
                double before = euler[2];
                euler[2] += deltatheta(3 * i + 2) * 180.0 / std::acos(-1) ;
                euler[2] = std::clamp(euler[2], double(b->rzmin), double(b->rzmax));
                if (euler[2] != before) limitExceeded = true; // bonus
            }
        }
    }
    // TODO#3-5
    // Return whether IK is stable
    // i.e. whether the ball is reachable
    // Hint:
    //      1. comment out the line here and return whether the IK is stable or not
    //      2. if the ball is reachable,  swinging its hand in air
    forwardSolver(posture, root_bone);
    double finalDist = (end_bone->end_position - target_pos).head<3>().norm();
    bool reachable = (finalDist < epsilon) && !limitExceeded;
    
    if (!reachable) {
        posture = original_posture;
        forwardSolver(posture, root_bone);
        return false;
    } 

}

}  // namespace kinematics
