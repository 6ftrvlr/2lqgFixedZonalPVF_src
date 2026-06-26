Original 2LQG:
$$\frac{\partial Q_i}{\partial t} + [\Psi_i, Q_i] = \nu_i\nabla^2Q_i$$

Decomposing into zonal mean and perturbation around it:
$$Q_i = \overline{Q}_i + q_i,\qquad\Psi_i = \overline{\Psi}_i + \psi_i$$

Zonally averaging the original 2LQG equation, using the formula $\overline{[A,B]}=\partial_y \overline{AB_x}$:
$$\frac{\partial\overline{Q}_i}{\partial t} = \frac{\partial}{\partial y}\left(\nu_i\,\partial_y\overline{Q}_i - \overline{(\Psi_i)_xQ_i}\right)\quad$$

considering
$$\partial_ y\overline{\Psi_x Q} = \partial_ y\overline{(\overline{\Psi}_x+\psi_x)(\overline{Q}+q)} = \partial_y\overline{\psi_x q}$$

therefore
$$\frac{\partial\overline{Q}}{\partial t} + \frac{\partial\overline{\psi_xq}}{\partial y} = \frac{\partial(\nu\partial_y\overline{Q})}{\partial y}$$

If the zonally average of $Q$ reaches a statistically steady state, then this gives an initial integral, which is 
$$\nu_1\,\partial_y\overline{Q_1} - \overline{\psi_1q_1} = -F_1,\qquad \nu_2\,\partial_y\overline{Q_2} - \overline{\psi_2q_2} = -F_2.$$

Again, from the original 2LQG equation, we derive the evolution of the perturbation field. Inserting the decomposition into the original equation (subscript $i$ omitted):
$$\frac{\partial(\overline{Q}+q)}{\partial t} + [\overline{\Psi}+\psi, \overline{Q}+q] = \nu\nabla^2(\overline{Q}+q)$$

i.e. (partially derivating a zonally averaged quantity wrt $x$ yields zero)
$$\frac{\partial\overline{Q}}{\partial t} + \frac{\partial q}{\partial t} + \psi_x(\overline{Q}_y+q_y) - (\overline{\Psi}_y+\psi_y)q_x = \nu_i\partial_{yy}\overline{Q} + \nu_i\nabla^2q$$

inserting the evolution of the zonally averaged PV, one derives
$$\frac{\partial q}{\partial t} + (\psi_x\overline{Q_y} - \overline{\Psi}_yq_x) + [\psi, q] - \frac{\partial\overline{\psi_x q}}{\partial y} = \nu\Delta q.$$

using again the formula $\overline{[A,B]}=\partial_y \overline{AB_x}$, this evolution equation is no more than
$$\frac{\partial q}{\partial t} + (\psi_x\overline{Q_y} - \overline{\Psi}_yq_x) + [\psi, q] - \overline{[\psi, q]} = \nu\Delta q.$$

Now one can assign a base state which fits in the original equation. In the BCI explorations, one usually selects
$$\psi_1 = -U y,\qquad\psi_2 = U y,\qquad Q_1 = k_D^2U y,\qquad Q_2 = -k_D^2U y.$$

Oh! I found where the problem is. Again, I assumed that $U=U(t)$. If not so, the base state fits the original equation. But now $U=U(t)$, the base state does not fit the equation. But in the traditions of the research into wave-mean-flow interaction, one may assume that $U'(t)$ is small enough so that one can neglect its influence. But is it ok this time? What needs to be further justified? In fact it is OK if we justify our work by verifying a quasi-static assumption concealed in the model shift.

We insert this background field and get our governing equation
$$\begin{aligned}
\frac{\partial q_1}{\partial t} &+ U[(q_1)_x+k_D^2(\psi_1)_x] + [\psi_1, q_1] - \overline{[\psi_1, q_1]} = \nu\Delta q_1,\\
\frac{\partial q_2}{\partial t} &- U[(q_2)_x+k_D^2(\psi_2)_x] + [\psi_2, q_2] - \overline{[\psi_2, q_2]} = \nu\Delta q_2,\\
F_1 &= \overline{(\psi_1)_xq_1} - \nu k_D^2 U,\\
F_2 &= \overline{(\psi_2)_xq_2} + \nu k_D^2 U.
\end{aligned}$$

One can now push a zonal flux constraint to 1st layer (or both?) by using the formulas and setting $F_1$ to be a constant, and therefore we constructs our problem.