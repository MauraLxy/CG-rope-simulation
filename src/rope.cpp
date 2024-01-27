#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {
    // Add a global variable damping_factor
    float damping_factor = 0.00005;

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Task 1.2): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.
        // Get the distance between each node.
        Vector2D distance = (end - start)/(num_nodes-1);
        // Set a mass with last position
        Mass* m1 = nullptr;
        // Iterate through all the nodes
        for (int index = 0; index < num_nodes; index++)
        {
            // Find the position of the next node
            Vector2D position = start + distance * index;
            // Create a new mass
            bool pinned = false;
            Mass* m2 = new Mass(position, node_mass, pinned);
            masses.push_back(m2);
            // Check whether the node is pinned
            for (auto &i : pinned_nodes) {
               masses[i]->pinned = true;
            }
            // Create a new spring
            if (m1 != nullptr){
                Spring* spring = new Spring(m1, m2, k);
                springs.push_back(spring);
            }
            // Update the previous position
            m1 = m2;
        }
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        bool is_explicit = false; // Switch between explicit and implicit Euler

        for (auto &s : springs)
        {
            // TODO (Task 1.3): Use Hooke's law to calculate the force on a node
            Vector2D direction = s->m2->position - s->m1->position;
            float length = direction.norm();
            Vector2D force = -s->k * direction/length * (length - s->rest_length);
            // add forces 
            s->m1->forces -= force;
            s->m2->forces += force;
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Task 1.6): Add global damping
                // use f_d = -k_d * v, where f_d is the damping force, k_d is a damping constant, and v is the velocity of the mass
                Vector2D damping = -damping_factor * m->velocity;
                m->forces += damping;
                // Add the gravity
                m->forces += gravity;
                // F = ma
                Vector2D acceleration = m->forces / m->mass;

                if (is_explicit) {
                    // TODO (Task 1.3): Add the force due to gravity, then compute the new velocity and position (explicit Euler)
                    // Formula: x(t+1) = x(t) + v(t) * dt
                    m->position += m->velocity * delta_t;
                    // Formula: v(t+1) = v(t) + a(t) * dt
                    m->velocity += acceleration * delta_t;
                } else {
                    // TODO (Task 1.3): Add the force due to gravity, then compute the new velocity and position (semi-implicit Euler)
                    // Formula: v(t+1) = v(t) + a(t) * dt
                    m->velocity += acceleration * delta_t;
                    // Formula: x(t+1) = x(t) + v(t+1) * dt
                    m->position += m->velocity * delta_t;
                }
            }
            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Task 1.5): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            Vector2D direction = s->m2->position - s->m1->position;
            float length = direction.norm();
            Vector2D force = -s->k * direction/length * (length - s->rest_length);
            // add forces 
            s->m1->forces -= force;
            s->m2->forces += force;
            // Correction vector sould be proportional to the displacement between the two masses
            float proportion = (length - s->rest_length) / length;
            // Correction vector should be in direction of one mass point to the other
            Vector2D correction = direction * proportion;

            // Update positions of masses to maintain the spring length.
            // Each mass should be moved by half of the displacement excluding pinned
            if(!s->m1->pinned && !s->m2->pinned){
                s->m1->position += 0.5f * correction;
                s->m2->position -= 0.5f * correction;
            } else if(!s->m1->pinned && s->m2->pinned) {
                s->m1->position += correction;
            } else if (!s->m2->pinned && s->m1->pinned){
                s->m2->position -= correction;
            }
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                Vector2D temp_position = m->position;
                // TODO (Part 1.5): Set the new position of the rope mass
                // Add gravity
                m->forces += gravity;
                // Formula: x(t+1) = x(t) + [x(t)-x(t-1)] + a(t) * dt * dt
                m->position += (temp_position - m->last_position) + m->forces/m->mass * delta_t * delta_t;

                // TODO (Part 1.6): Add global Verlet damping
                // Formula: x(t+1) = x(t) + (1-damping_factor) * [x(t)-x(t-1)] + a(t) * dt * dt
                // m->position += (1 - damping_factor) * (m->position-m->last_position) + m->forces/m->mass * delta_t * delta_t;

                // Update the previous position to be the current position
                m->last_position = temp_position; 
            }
            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }
}
