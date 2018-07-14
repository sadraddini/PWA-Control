"""
Created on Mon June 25 2018

@author:    Sadra Sadraddini
            CSAIL, MIT 
            32 Vassar st, Cambridge, MA 02139
            s a d r a @ mit.edu
"""


# Primary imports
import numpy as np
from gurobipy import Model,GRB,LinExpr,QuadExpr
from random import choice as rchoice

# Secondary imports
from auxilary_methods import find_mode,valuation,mode_sequence
from ana_system import state

def polytope_trajectory(s,x0,state_end,T,alpha_start):
    model=Model("trajectory of polytopes")
    x={}
    u={}
    theta={}
    z={}
    G_bound=100
    # Mode 1:
    for t in range(T):
        x[t]=np.empty((s.n,1),dtype='object') # n*1
        u[t]=np.empty((s.m,1),dtype='object') # m*1
        theta[t]=np.empty((s.m,s.n),dtype='object') # n*m
        for row in range(s.n):
            x[t][row,0]=model.addVar(lb=-G_bound,ub=G_bound)
        for row in range(s.m):
            u[t][row,0]=model.addVar(lb=-G_bound,ub=G_bound)
        for row in range(s.m):
            for column in range(s.n):
                theta[t][row,column]=model.addVar(lb=-G_bound,ub=G_bound)   
    for t in range(T+1):
        for i in s.modes:
            z[t,i]=model.addVar(vtype=GRB.BINARY)
    x[T]=np.empty((s.n,1),dtype='object') # Final state in Mode i
    for row in range(s.n):
        x[T][row,0]=model.addVar(lb=-G_bound,ub=G_bound)
    G={}
    for t in range(T+1):
        G[t]=np.empty((s.n,s.n),dtype='object')
        for row in range(s.n):
            for column in range(s.n):
                G[t][row,column]=model.addVar(lb=-G_bound,ub=G_bound)
    model.update()
    # Trajectory Constraints:
    # Mode i:
    bigM=100
    for i in s.modes:
        for t in range(T):
            for row in range(s.n):
                Ax=LinExpr()
                for k in range(s.n):
                    Ax.add(s.A[i][row,k]*x[t][k,0])
                for k in range(s.m):
                    Ax.add(s.B[i][row,k]*u[t][k,0])
                model.addConstr(x[t+1][row,0]<=Ax+s.c[i][row]+bigM-bigM*z[t,i])
                model.addConstr(x[t+1][row,0]>=Ax+s.c[i][row]-bigM+bigM*z[t,i])
    # Generator Dynamics Constraints:
    for i in s.modes:
        for t in range(T):
            for row in range(s.n):
                for column in range(s.n):
                    AG=LinExpr()
                    for k in range(s.n):
                        AG.add(s.A[i][row,k]*G[t][k,column])
                    for k in range(s.m):
                        AG.add(s.B[i][row,k]*theta[t][k,column])
                    model.addConstr(G[t+1][row,column]<=AG+bigM-bigM*z[t,i])
                    model.addConstr(G[t+1][row,column]>=AG-bigM+bigM*z[t,i])
    # Constraints of modes:
    for t in range(T+1):
        sum_z=LinExpr()
        for i in s.modes:
            sum_z.add(z[t,i])
        model.addConstr(sum_z==1)
    # Constraints of mode subsets
    for t in range(T):
        for i in s.modes:
            subset_MILP(model,G[t],s.Pi,s.H[i],s.h[i],x[t],z[t,i])
            subset_MILP(model,theta[t],s.Pi,s.F[i],s.f[i],u[t],z[t,i])   
    # set objective
    J_area=LinExpr()
    d_min=model.addVar(lb=0.0001)
    beta=10**1 # Weight of infinity norm
    model.update()
    for row in range(s.n):
        for column in range(s.n):
            if row<column:
#                continue
                model.addConstr(G[0][row,column]==0)
            elif row==column:
                model.addConstr(G[0][row,column]>=d_min)
    J_area.add(-d_min*T*s.n*beta)
    for row in range(s.n):
        for t in range(T):
            J_area.add(-G[t][row,row])
    # Terminal Constraint
    terminal_constraint(s,x,G,T,model,state_end)
    # Starting Point
    i_start=find_mode(s,x0)
    for i in s.modes:
        model.addConstr(z[0,i]==int(i==i_start))
#    model.setParam('OutputFlag',False)
    if alpha_start==-1:
        eps=0.1
        x_delta={}
        for row in range(s.n):
            x_delta[row]=model.addVar(lb=-eps,ub=eps)
        model.update()
        print("only area!")
        for row in range(s.n):
            model.addConstr(x[0][row,0]==x0[row,0]+x_delta[row])
        model.setObjective(J_area)
        model.optimize()
    else:          
        J_start=QuadExpr()
        for row in range(s.n):
            J_start.add(s.weight[row]*x[0][row,0]*x[0][row,0]-2*s.weight[row]*x0[row,0]*x[0][row,0])
        model.setObjective(J_area+alpha_start*J_start)
        model.optimize()
    if model.Status!=2 and model.Status!=11:
        flag=False
        print("*"*20,"Flag is False and Status is",model.Status)
        return (x,u,G,theta,z,flag)
    else:
        flag=True
#        print("*"*20,"Flag is True and Status is",model.Status)
        x_n=valuation(x)
        u_n=valuation(u)
        G_n=valuation(G)
        theta_n=valuation(theta)
        z_n=mode_sequence(s,z)
#        print("--"*20,"final")
#        print ("new:",G_n[0],np.linalg.det(G_n[0]))
        print("initial x0 was",x0.T)
        print("final x0 is",x_n[0].T)
        print("final G0 is",G_n[0])
        if abs(np.linalg.det(G_n[0]))<10**-5:
            flag=False
        return (x_n,u_n,G_n,theta_n,z_n,flag)
    
def grow_tree_exact(s,i,T,alpha,eps=0.01):
    x0=sample(s.l[i],s.u[i])
    if is_inside_machine_state(s,x0,eps)==True:
        print("covered!: %d"%i)
        return
    else:
        for branch in list(s.branches)[::-1]:
            for state_x in list(branch):
#                print("evaluating with alpha",state_x,"alpha=",alpha)
                for Tau in range(5,T):
                    s.factor=1.0
                    (x,u,G,theta,z,flag)=funnel_to_state_multi_mode(s,x0,state_x,Tau,alpha)
                    x0_new=x[0]
                    if flag==True:
                        if are_all_vertices_inside_machine(s,x0_new,G[0],eps)==True:
                            print("all vertices covered after optimization! tau=",Tau,"x0 was",x0.T," x0_new=",x0_new.T)
                            continue
                        else:
                            make_state_trajectory_state_end(s,x,u,z,G,theta,Tau,state_x)
                            return # Ignorring the sript below:
                        if is_inside_machine_state(s,x0_new,eps)==True:
                            print("covered after optimization! tau=",Tau,"x0 was",x0.T," x0_new=",x0_new.T)
                            continue
                        else:
                            print("generating exactness",state_x)
                            (x,u,G,theta,z,flag)=funnel_to_state_multi_mode(s,x0_new,state_x,Tau,10000)
                            if flag==False:
                                print("flag is false! for Tau=",Tau)
                                continue
                                raise(ValueError("flag is false! for Tau=",Tau))
                            make_state_trajectory_state_end(s,x,u,z,G,theta,Tau,state_x)
                            return
                    else:
                        print("no %d trajectory exists to state with"%Tau ,state_x.x.T,"alpha was",alpha)
                        continue    

def grow_tree_exact_initial(s,x0,T,alpha,eps=0.01,max_iterations=10):
    for branch in list(s.branches)[::-1]:
        for state_x in list(branch):
            (x,u,G,theta,z,flag)=funnel_to_state_multi_mode(s,x0,state_x,T,alpha)
            if flag==True:
                print("flag is TRUE!")
                if are_all_vertices_inside_machine(s,x[0],G[0],eps)==True:
                    print("all vertices covered after optimization! T=",T,"x0 was",x0.T," x0_new=",x[0].T)
                    continue
                else:
                    print("goal state is",state_x)
                    make_state_trajectory_state_end(s,x,u,z,G,theta,T,state_x)
                    return 
            else:
                print("no %d trajectory exists to state with"%T ,state_x.x.T,"alpha was",alpha)
                continue
    return                  
    
def are_all_vertices_inside_machine(s,x,G,eps):
    v=vertices_cube(G.shape[1])
    vertices=(np.dot(G,v.T)+x).T
    if  is_inside_machine_state(s,x,eps)==False:
        return False
    for i in range(vertices.shape[0]):
        vertex_to_check=vertices[i,:].reshape(G.shape[1],1)
        if is_inside_machine_state(s,vertex_to_check,eps)==False:
            return False
    return True
                
def make_state_trajectory_free(s,x,u,seq,G,theta,T):
    x_temp={}
    ID=rchoice(range(1000))
    for t in range(0,T+1):
        x_temp[t]=state(x[t],G[t],seq[t],ID,t,6)
    # Build Transitons
    for t in range(0,T):
        x_temp[t].edge.append((x_temp[t+1],u[t],theta[t]))
    s.X.extend(x_temp.values())
    s.branches.append(x_temp.values())
    
def make_state_trajectory_state_end(s,x,u,seq,G,theta,T,state_end):
    x_temp={}
    ID=rchoice(range(1000))
    for t in range(0,T):
        x_temp[t]=state(x[t],G[t],seq[t],ID,t,7)
        x_temp[t].time_to_go=T-t+state_end.time_to_go
    # Build Transitons
    for t in range(0,T-1):
        x_temp[t].successor=(x_temp[t+1],u[t],theta[t])
    x_temp[T-1].successor=(state_end,u[T-1],theta[T-1])
    s.X.extend(x_temp.values())
    

    
def make_zero_alpha_backward_state(s,state_goal,T,i,eps):
    x0=sample(s.l[0],s.u[0])
    (x,u,G,theta,z,flag)=funnel_to_state_multi_mode(s,x0,state_goal,T,0)
    if flag==True:
        print("flag is TRUE!")
        if are_all_vertices_inside_machine(s,x[0],G[0],eps)==True:
            print("all vertices covered after optimization! T=",T,"x0 was",x0.T," x0_new=",x[0].T)
            state_goal.backward_zero=-1
        else:
            make_state_trajectory_state_end(s,x,u,z,G,theta,T,state_goal)
            state_goal.backward_zero=1
            
    else:
        print("no alpha-zero %d trajectory exists to state with"%T ,state_goal.x.T)
        state_goal.backward_zero=-2

def make_zero_alpha_backward_set(s,T,i,eps):
    for state_goal in [x for x in s.X if x.mode==i and x.mode==0]:
        make_zero_alpha_backward_state(s,state_goal,T,i,eps)

def subset_MILP(model,G,Pi,H,h,x,z_time_mode):
    """
    Description: Add Farkas lemma constraints for subset inclusion of x+GP in {e|H.<h}
    Inputs: 
        model: Gurobi optimization model
        G: n * n_g generator matrix
        Pi:  n_pi*n matrix where {x|Pi x< 1} is the primitive polytope
        {x| Hx<h} is the set constraint
        x is the point
    Output:
        no direct output. Adds constraints to the model. 
        FUTURE: we may add lambda positive matrix as an output for debugging
    """
    bigM=100
    (n,n_g)=G.shape
    (n_p,n)=Pi.shape
    (n_h,n)=H.shape
    Lambda=np.empty((n_h,n_p),dtype='object')
    for row in range(n_h):
        for column in range(n_p):
            Lambda[row,column]=model.addVar(lb=0)
    model.update()
    # Lambda * Pi = H * G
    for row in range(n_h):
        for column in range(n_g):
            s_left=LinExpr()
            s_right=LinExpr()
            for k in range(n_p):
                s_left.add(Lambda[row,k]*Pi[k,column])
            for k in range(n):
                s_right.add(H[row,k]*G[k,column])
            model.addConstr(s_left==s_right)
    # Lambda * 1 <= H*
    for row in range(n_h):
        s_left=LinExpr()
        s_right=LinExpr()
        for k in range(n_p):
            s_left.add(Lambda[row,k])
        for k in range(n):
            s_right.add(H[row,k]*x[k,0])
        model.addConstr(s_left<=h[row,0]-s_right+bigM-bigM*z_time_mode) 

def subset(model,G,Pi,H,h,x):
    """
    Description: Add Farkas lemma constraints for subset inclusion of x+GP in {e|H.<h}
    Inputs: 
        model: Gurobi optimization model
        G: n * n_g generator matrix
        Pi:  n_pi*n matrix where {x|Pi x< 1} is the primitive polytope
        {x| Hx<h} is the set constraint
        x is the point
    Output:
        no direct output. Adds constraints to the model. 
        FUTURE: we may add lambda positive matrix as an output for debugging
    """
    (n,n_g)=G.shape
    (n_p,n)=Pi.shape
    (n_h,n)=H.shape
    Lambda=np.empty((n_h,n_p),dtype='object')
    for row in range(n_h):
        for column in range(n_p):
            Lambda[row,column]=model.addVar(lb=0)
    model.update()
    # Lambda * Pi = H * G
    for row in range(n_h):
        for column in range(n_g):
            s_left=LinExpr()
            s_right=LinExpr()
            for k in range(n_p):
                s_left.add(Lambda[row,k]*Pi[k,column])
            for k in range(n):
                s_right.add(H[row,k]*G[k,column])
            model.addConstr(s_left==s_right)
    # Lambda * 1 <= H*
    for row in range(n_h):
        s_left=LinExpr()
        s_right=LinExpr()
        for k in range(n_p):
            s_left.add(Lambda[row,k])
        for k in range(n):
            s_right.add(H[row,k]*x[k,0])
        model.addConstr(s_left<=h[row,0]-s_right) 

def terminal_constraint(s,x,G,T,model,state):
        # Terminal Constraint
    if state.character!=-1 and state.volume_flag==True:
        print("taking G approach")
        H=np.dot(s.Pi,state.Ginv)
        h=np.ones((s.Pi.shape[0],1))+np.dot(H,state.x)
        #print("Nonzero Volume!","Ginv=",state2.Ginv,"H=",H,"h=",h)
        subset(model,G[T],s.Pi,H,h,x[T])
    else:
        print("taking Lambda approach")
        # The vertices by G[i,T]
        Lambda={}
        for alpha in range(2**s.n):
            for beta in range(2**state.G.shape[1]):
                Lambda[alpha,beta]=model.addVar(lb=0)
        model.update()
        for alpha in range(2**s.n):
            for row in range(s.n):
                exp_left=LinExpr()
                exp_right=LinExpr()
                for k in range(s.n):
                    exp_left.add(G[T][row,k]*s.vertices[alpha,:].reshape(s.n,1)[k,0])
                for beta in range(2**state.G.shape[1]):
                    exp_right.add(Lambda[alpha,beta]*state.vertices[beta,:].reshape(s.n,1)[row,0])
                model.addConstr(x[T][row,0]+exp_left==exp_right+state.x[row,0])
            lambda_sum=LinExpr()
            for beta in range(2**state.G.shape[1]):
                lambda_sum.add(Lambda[alpha,beta])
            model.addConstr(lambda_sum<=1)