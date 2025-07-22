from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from openmm import XmlSerializer
from parmed.amber import LoadParm
from parmed.openmm import RestartReporter
from parmed import unit as u
from tqdm import tqdm


#Choosing the platform to run on
platform = Platform.getPlatformByName('CUDA')  
TOTAL_REPS = 16
for rep in tqdm(range(TOTAL_REPS)):
    parm7_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/parm7/MAAM_{rep}.parm7'
    input_rst_path =f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/equil/MAAM_{rep}.rst7'
    dcd_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/prod/MAAM_{rep}.dcd'
    log_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/prod/MAAM_{rep}.log'
    chk_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/prod/MAAM_{rep}.chk'
    rst_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/prod/MAAM_{rep}.rst7'
    #load_chk_path = f'/home/hassan/PycharmProjects/peptide_conduction/tetra_peptides/npt/chk/{i}_NPT_final.chk'
    equil_state_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/prod/MAAM_{rep}_final.state'
    equil_chk_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/prod/MAAM_{rep}_final.chk'

    #print(f"Running production run for 10ns with e-field for {i} system")
    
    parm=LoadParm(parm7_path, input_rst_path)
    mdsteps = 100000000
    
    dis = 0.6

    #System Parameters
    nonbondedMethod = PME
    nonbondedCutoff = 1.2*nanometers
    ewaldErrorTolerance = 0.0005
    constraints = HBonds
    rigidWater = True
    constraintTolerance = 0.000001

    #Simulation parameters
    dt = 0.002*picoseconds
    temperature = 300*kelvin
    friction = 2.8284/picosecond
    pressure = 1.0*atmospheres


    dcdReporter = DCDReporter(dcd_path, 50000)
    dataReporter = StateDataReporter(log_path, 50000, totalSteps= mdsteps,
        step=True, speed=True, progress=True, elapsedTime=True, remainingTime=True, potentialEnergy=True, temperature=True, volume=True, density=True, separator=',')
    checkpointReporter = CheckpointReporter(chk_path, 50000)
    restartreporter = RestartReporter(rst_path,reportInterval=50000,netcdf=True)
    


    topology = parm.topology
    positions = parm.positions
    n_atoms = len(positions)
    system = parm.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
    integrator = LangevinIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)
    #========================================================================#

    s_inds = []
    cg_inds = []
    ce_inds = []
    for i, _ in enumerate(parm.residues):
        if parm.residues[i].name =='MET':
            for j in enumerate (parm.residues[i].atoms):
                if j[1].name == 'SD':
                    s_inds.append(j[1].idx)
                if j[1].name == 'CG':
                    cg_inds.append(j[1].idx)
                if j[1].name == 'CE':
                    ce_inds.append(j[1].idx)
                
    s_ind1 = s_inds[0]
    s_ind2 = s_inds[1]
    
    cg_ind1 = cg_inds[0]
    cg_ind2 = cg_inds[1]
    
    ce_ind1 = ce_inds[0]
    ce_ind2 = ce_inds[1]



    #================================================================================#
    # add S-S force via CustomCompoundBondForce
    sspull_constant    = 418.4 # units: kJ/(mol*nm^2)
    sspull_target_dist = dis   # units: nm
    sspull = CustomCompoundBondForce(2, '0.5*k*((z2 - z1) - targetdist)^2')
    sspull.addGlobalParameter('k', sspull_constant)
    sspull.addGlobalParameter('targetdist', sspull_target_dist)
    sspull.addBond([s_ind1, s_ind2], [])
    system.addForce(sspull)

   #========== Force orienting sulfur lone pair to z- axis ===========================# 


    ssorient_constant = 4.184 * 10. # 10 kcal/mol when unaligned, -10 kcal/mol when aligned
    ssorient = CustomCompoundBondForce(3, 'dir*k*(z1 - (z2+z3)/2)/pointdistance(x1, y1, z1, (x2+x3)/2, (y2+y3)/2, (z2+z3)/2)')
    ssorient.addPerBondParameter('dir')
    ssorient.addGlobalParameter('k', ssorient_constant)
    ssorient.addBond([s_ind1, cg_ind1, ce_ind1], [1.])
    ssorient.addBond([s_ind2, cg_ind2, ce_ind2], [-1.])
    system.addForce(ssorient)

 
    #==========================================================#
    # add electric field
    applied_potential   = 0.00025       # Volts
    separation_distance = sspull_target_dist  # nm
    au_s_bond_length    = 0.25      # nm
    avogadro_number     = 6.022e23
    electron_volt       = 1.602e-19 # C/e-
    J_to_kJ             = 0.001
    efield_strength     = applied_potential * electron_volt * avogadro_number * J_to_kJ / (separation_distance + 2. * au_s_bond_length)
    efield = CustomExternalForce('q*Ez*z')
    efield.addPerParticleParameter('q')
    efield.addGlobalParameter('Ez', efield_strength) # kJ/(mol*nm*e-)
    nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
    for i in range(n_atoms):
        charge, sigma, epsilon = nonbonded.getParticleParameters(i)
        efield.addParticle(i, [charge])
    system.addForce(efield)


    #=========================================================#

    simulation = Simulation(topology, system, integrator, platform)#, platformProperties)
    if parm.box_vectors is not None:
        simulation.context.setPeriodicBoxVectors(*parm.box_vectors)


    #simulation.loadCheckpoint(load_chk_path)
    #NPT_state = simulation.context.getState(getVelocities=True, getPositions=True)
    #positions = NPT_state.getPositions()
    #velocities = NPT_state.getVelocities()
    positions=parm.positions
    simulation.context.setPositions(positions)
    simulation.context.setVelocitiesToTemperature(temperature)

    simulation.reporters.append(dcdReporter)
    simulation.reporters.append(dataReporter)
    simulation.reporters.append(checkpointReporter)
    simulation.reporters.append(restartreporter)

    barostat = system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin))
    simulation.context.reinitialize(True)
    simulation.currentStep = 0
    simulation.step(mdsteps)
    simulation.saveState(equil_state_path)
    simulation.saveCheckpoint(equil_chk_path)



