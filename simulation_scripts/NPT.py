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
names = ["MAAM"]
TOTAL_REPS = 16
for rep in tqdm(range(TOTAL_REPS)):
    parm7_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/parm7/MAAM_{rep}.parm7'
    crd_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/rst7/MAAM_{rep}.rst7' 
    dcd_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/npt/MAAM_{rep}.dcd'
    log_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/npt/MAAM_{rep}.log'
    chk_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/npt/MAAM_{rep}.chk'
    rst_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/npt/MAAM_{rep}_restart.rst7'
    input_rst_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/min/MAAM_{rep}_restart.rst7'
    load_chk_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/min/MAAM_{rep}.chk'
    npt_state_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/npt/MAAM_{rep}_final.state'
    npt_chk_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/npt/MAAM_{rep}_final.chk'

    
    parm=LoadParm(parm7_path, input_rst_path)
    mdsteps = 250000
    

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
    dataReporter = StateDataReporter(log_path, 50000, totalSteps=2*mdsteps,
        step=True, speed=True, progress=True, elapsedTime=True, remainingTime=True, potentialEnergy=True, temperature=True, volume=True, density=True, separator=',')
    checkpointReporter = CheckpointReporter(chk_path, 50000)
    restartreporter = RestartReporter(rst_path,reportInterval=50000,netcdf=True)
    


    topology = parm.topology
    positions = parm.positions
    system = parm.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
    integrator = LangevinIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)


    totalrestrainedatoms=['N','CA','C']
    force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k", 1.0*kilocalories_per_mole/angstroms**2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    c=0
    for i, atom_crd in enumerate(parm.positions):
        if parm.atoms[i].name in totalrestrainedatoms:
            c+=1
            force.addParticle(i, atom_crd.value_in_unit(u.nanometers))
    system.addForce(force)


    simulation = Simulation(topology, system, integrator, platform)#, platformProperties)
    if parm.box_vectors is not None:
        simulation.context.setPeriodicBoxVectors(*parm.box_vectors)


    simulation.loadCheckpoint(load_chk_path)
    NVT_state = simulation.context.getState(getVelocities=True, getPositions=True)
    positions = NVT_state.getPositions()
    velocities = NVT_state.getVelocities()
    simulation.context.setPositions(positions)
    simulation.context.setVelocities(velocities)

    simulation.reporters.append(dcdReporter)
    simulation.reporters.append(dataReporter)
    simulation.reporters.append(checkpointReporter)
    simulation.reporters.append(restartreporter)
	
	
	
    #print(f"Running NPT for 2ns for {j} system")
    barostat = system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin))
    simulation.context.reinitialize(True)
    simulation.currentStep = 0
    simulation.step(mdsteps)
    #print('Running NPT equilibration...')
    for i in range(100):
      simulation.step(int(mdsteps/100))
      #print(f'k = {float(99.02/100-(i*0.98/100))}')
      simulation.context.setParameter('k', (float(99.02/100-(i*0.98/100))*kilocalories_per_mole/angstroms**2))
    simulation.context.setParameter('k', 0)
    simulation.step(mdsteps)
    # simulation.minimizeEnergy()
    # simulation.saveState('minimized.state')
    # simulation.saveCheckpoint('minimized.chk')

    # simulation.context.setVelocitiesToTemperature(5*kelvin)
    # print('Warming up the system...')
    # T = 5
    # for i in range(60):
    #   simulation.step(int(mdsteps/60))
    #   temperature = (T+(i*T))*kelvin
    #   integrator.setTemperature(temperature)

    simulation.saveState(npt_state_path)
    simulation.saveCheckpoint(npt_chk_path)



