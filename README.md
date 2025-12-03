# Interstellar

I am building a multi-satellite swarm system based on principles of Flocking and Cyclic Pursuit. This work is a subset of a bigger project titled **_Swarm-Based Encirclement of Orbital Debris in Low Earth Orbit_**. You can find a preliminary 'white paper' like guide to the project at the bottom.
Why 'Interstellar'? The notion of ambitious space missions piloted by humans is a real and not just cinematic marvel. Christopher Nolan does a great job of portraying this of Cooper in his 2014 Oscar Nominee - Interstellar. However, it always got me wondering - what if TARS and CASE took over the mission entirely? So here's an extension of what graduate studies in Space Robotics and a deep admiration for cool space tech led to - start with space debris (bigger problem that fewer people care about), and solve autonomy incrementally. 

### Mission Overview

Fragmentation events relating to collisions, fuel, and anti-satellite tests and more, have contributed
to a significant increase in actively tracked objects in Low Earth Orbit. Since 2024, over 3000 more
tracked objects were introduced due to fragmentation1. That said, we propose to design and simulate
a cooperative satellite swarm system composed of five autonomous nano-satellites operating in Low
Earth Orbit which will coordinate to track and encircle orbital debris in a cage-like, circular formation.
They will do so while maintaining their own orbital paths around Earth.

The swarm’s relative motion and coordination will be based on Olfati-Saber Flocking (OSF) and
cyclic pursuit principles to ensure collision avoidance, cohesion, and stable encirclement of the target.
The system will be implemented and tested in MuJoCo, serving as a physics-based environment to
validate dynamics, tracking performance, and the feasibility of containment strategies for active debris
mitigation.

### Experiments
I stage development by testing algorithmic improvements in MATLAB, and then in MuJoCo. Below is a table of cyclic pursuit algorithms I have tested and adapted for this project. The algorithms operate in Earth-Centred Inertial (ECI) frame of reference - good for orbital mechanics i.e., HCW and pursuit laws assume a non-rotating inertial frame.

| File Name                       | Concept Used                               | 2D/3D | Linear / Non-linear | Frame of Reference | Final Control Law \(u_i\) (compact form)                                                               |
|---------------------------------|--------------------------------------------|-------|---------------------|--------------------|--------------------------------------------------------------------------------------------------------|
| cyclic_pursuit_velocity_field.m | Single-integrator velocity field + pursuit | 2D    | Linear              | Inertial-like      | $u_i = -k_v(v_i - v_{\text{des},i})\) where \(v_{\text{des},i} = k_d R(\alpha)(x_{i+1}-x_i)-k_c x_i$   |
| cyclic_pursuit_direct_acc.m     | Direct double-integrator pursuit           | 2D    | Linear              | Inertial-like      | $u_i = k_d R(\alpha)(x_{i+1}-x_i) + R(\alpha)(v_{i+1}-v_i) - (k_ck_d)x_i - (k_c+k_d)v_i$               |
| cyclic_pursuit_velocity3D.m     | 3D velocity-field pursuit                  | 3D    | Linear              | Inertial-like      | $u_i = -k_v\,(v_i - v_{\text{des},i})\), \(v_{\text{des},i} = k_d R_z(\alpha)(r_{i+1}-r_i) - k_c r_i$  |
| cyclic_pursuit_adapt.m          | CW-inspired pursuit (no radius keeping)    | 3D    | Linearized CW       | LVLH-like          | $u_i = -f_i + k_g( k_dTRT^{-1}(r_{i+1}-r_i) + TRT^{-1}(v_{i+1}-v_i) - k_ck_d r_i - (k_c+k_d/k_g)v_i )$ |
| cyclic_pursuit_adapt_radius.m   | CW pursuit + radial PD                     | 3D    | Linearized CW       | LVLH-like          | $u_i = -f_i + k_g(\text{pursuit terms}) + u_r\), where \(u_r = -k_r e_r n_i - k_{vr} v_r n_i$          |
| cyclic_pursuit_cw.m             | Full CW-style pursuit (pos/vel shaping)    | 3D    | Linearized CW       | LVLH-like          | $u_i = -f_i + k_g(k_dTRT^{-1}(r_{i+1}-r_i) + TRT^{-1}(v_{i+1}-v_i) + \text{spr} + \text{damp})$        |
> Note: CW stands for [Clohessy-Wiltshire](https://ensatellite.com/hills-equations/). It's a linearized system that captures the motion of a rigid body in an orbital frame. 

### Recommendations for Use

- **cyclic_pursuit_velocity_field.m**:
  Use when you want simple **2D cyclic pursuit** with clean circular motion.
  Good for **teaching**, **debugging**, and **visual intuition**.
  Best when you do *not* care about dynamics—just directional fields.

- **cyclic_pursuit_direct_acc.m**:
  Use when you need **double-integrator physics** (acceleration-level control).
  Good for UAV/robot simulations where agents have real inertia.
  Best for *non-orbital*, flat-ground experiments.

- **cyclic_pursuit_velocity3D.m**: 
  Use when you want **3D cyclic pursuit** without orbital mechanics.
  Good for drone swarms, multi-agent 3D formation tests.
  Use if you just need geometry, not orbital physics.

- **cyclic_pursuit_adapt.m**: 
  Use when exploring **Clohessy–Wiltshire-like interactions** but without radius enforcement.
  Suitable for early testing of CW-inspired pursuit ideas.
  Not recommended for stable long simulations (agents collapse inward).

- **cyclic_pursuit_adapt_radius.m**: 
  Use when you want **CW-style pursuit** *with* **stable radius keeping**.
  Best for simulating **relative motion near a nominal orbit**.
  Good precursor to integrating real orbital gravity in MuJoCo or ECI.

- **cyclic_pursuit_cw.m**: 
  Use when you want the closest behavior to **true orbital formation flying**.
  Best for **spacecraft swarm algorithms**, **formation keeping**, and **relative dynamics** research.
  Ideal if you plan to later add **ECI gravity + LVLH transformation**.

### Citation
I've largely used the following paper to draw theoretical validation from.

- Ramirez, J. L., et al. *Distributed Control of Spacecraft Formation via Cyclic Pursuit: Theory and Experiments.* American Control Conference (ACC), 2009, pp. 4811–4817. © 2009 IEEE.