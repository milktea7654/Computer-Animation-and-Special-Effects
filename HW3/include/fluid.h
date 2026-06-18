#pragma once

#include <glad/gl.h>
#include <vector>
#include <Eigen/Dense>

enum CellType
{
	AIR,
	FLUID,
	SOLID
};

class Fluid
{
public:
	Fluid(int num_particles, float radius, float obstacle_radius, float flip_ratio, float cell_dim, float relaxation_cell_dim, int iterations, float viewport_w, float viewport_h, float dt, float gravity, float stiffness, bool density_correction, float over_relaxation);
	~Fluid() = default;

	void update(int render_option, bool lmbDown);

	void setupRendering(GLuint pointShaders, GLuint obstacleShaders);
	void setObstacle(float mouse_x, float mouse_y, float mouse_vx, float mouse_vy);
	void setFlipRatio(float ratio) { flip_ratio = ratio; };
	void setStiffnessCoeff(float sc) { stiffness_coefficient = sc;  }
	void setDt(float dt_) { dt = dt_; };
	void setGravity(float grav) { gravity = grav; };
	void setDensityCorrection(bool dc) { density_correction = dc; };
	void setOverRelaxation(float or_) { over_relaxation = or_; };
	void setIterations(int iter) { iterations = iter; }
private:
	
	// General
	int iterations;
	float dt, gravity, viewport_w, viewport_h;
	float flip_ratio, stiffness_coefficient, over_relaxation;
	bool density_correction;

	// Particles
	std::vector<Eigen::Vector2f> particle_pos;
	std::vector<Eigen::Vector2f> particle_vel;
	float particle_radius;
	int num_particles;

	float particle_rest_density;

	// For relaxation only
	float relaxation_cell_dim;
	int relaxation_cell_rows, relaxation_cell_cols;
	std::vector<std::vector<std::vector<int>>> relaxation_cell_particle_ids;

	// Cells
	int cell_rows, cell_cols;
	float cell_dim;
	std::vector<std::vector<CellType>> cell_types;
	std::vector<std::vector<Eigen::Vector2f>> cell_velocities;
	std::vector<std::vector<Eigen::Vector2f>> prev_cell_velocities;
	std::vector<std::vector<float>> cell_densities;

	// obstacle
	float obstacle_x, obstacle_y, obstacle_vx, obstacle_vy, obstacle_r;

	// graphics
	GLuint point_shaders, obstacle_shaders,
		particles_vao, particles_vbo, particle_colors_vbo;
	GLuint cells_vao, cells_vbo, cell_colors_vbo;
	GLuint obstacle_vao;
	std::vector<Eigen::Vector2f> cell_centers_rendering;

	std::vector<Eigen::Vector3f> particle_colors;
	std::vector<Eigen::Vector3f> cell_colors;

	// physics
	void integrate();
	void particleRelaxation();
	void handleCollisions(bool lmbDown);
	void updateDensity();
	void transferVelocities(bool to_cell);
	void solveIncompressibility();

	// graphics
	void updateCellColors();
	void updateParticleColors();
	void renderParticles();
	void renderCells();
	void updateParticleBuffers();
	void updateCellColorBuffers();
	void renderObstacle();
};