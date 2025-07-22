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
    input_rst_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/npt/MAAM_{rep}_restart.rst7'
    dcd_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/equil/MAAM_{rep}.dcd'
    log_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/equil/MAAM_{rep}.log'
    chk_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/equil/MAAM_{rep}.chk'
    rst_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/equil/MAAM_{rep}.rst7'
    load_chk_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/npt/MAAM_{rep}.chk'
    equil_state_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/equil/MAAM_{rep}_final.state'
    equil_chk_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/equil/MAAM_{rep}_final.chk'

    #print(f"Running equilibration for 4ns for {j} system")
    
    parm=LoadParm(parm7_path, input_rst_path)
    mdsteps = 500000
    

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
    dataReporter = StateDataReporter(log_path, 50000, totalSteps=mdsteps,
        step=True, speed=True, progress=True, elapsedTime=True, remainingTime=True, potentialEnergy=True, temperature=True, volume=True, density=True, separator=',')
    checkpointReporter = CheckpointReporter(chk_path, 50000)
    restartreporter = RestartReporter(rst_path,reportInterval=50000,netcdf=True)
    


    topology = parm.topology
    positions = parm.positions
    system = parm.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
    integrator = LangevinIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)

    simulation = Simulation(topology, system, integrator, platform)#, platformProperties)
    if parm.box_vectors is not None:
        simulation.context.setPeriodicBoxVectors(*parm.box_vectors)


    simulation.loadCheckpoint(load_chk_path)
    NPT_state = simulation.context.getState(getVelocities=True, getPositions=True)
    positions = NPT_state.getPositions()
    velocities = NPT_state.getVelocities()
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



