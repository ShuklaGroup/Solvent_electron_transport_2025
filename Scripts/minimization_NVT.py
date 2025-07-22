import sys
sys.path.insert(0,'/home/hassan/research/anaconda3/envs/openmm2/lib/python3.10/site-packages')
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

names = ["MAAM"]
platform = Platform.getPlatformByName('CUDA')
TOTAL_REPS = 16
for rep in range(TOTAL_REPS):
    parm7_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/parm7/MAAM_{rep}.parm7'
    crd_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/rst7/MAAM_{rep}.rst7' 
    dcd_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/min/MAAM_{rep}.dcd'
    log_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/min/MAAM_{rep}.log'
    chk_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/min/MAAM_{rep}.chk'
    rst_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/min/MAAM_{rep}_restart.rst7'
    min_state_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/min/MAAM_{rep}.state'
    min_chk_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/min/MAAM_{rep}_minimized.chk'
    nvt_state_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/min/MAAM_{rep}_final.state'
    nvt_chk_path = f'/home/hassan/research/projects/peptide_environment/solvent_variability/acn/MAAM/sim/min/MAAM_{rep}_final.chk'

    parm=LoadParm(parm7_path , crd_path)

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
    #surfacetension=200*bar*nanometer
    #barostatInterval = 25
    mdsteps = 500000

    #equilibrationSteps = 100000
    


    
    dcdReporter = DCDReporter(dcd_path, 50000)
    dataReporter = StateDataReporter(log_path, 50000, totalSteps=mdsteps,
        step=True, speed=True, progress=True, elapsedTime=True, remainingTime=True, potentialEnergy=True, temperature=True, volume=True, density=True, separator=',')
    checkpointReporter = CheckpointReporter(chk_path, 50000)
    restartreporter = RestartReporter(rst_path,reportInterval=50000,netcdf=True)



    topology = parm.topology
    positions = parm.positions
    
    
    system = parm.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
    #system.addForce(MonteCarloMembraneBarostat(pressure,surfacetension,temperature,MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree))

    integrator = LangevinIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(constraintTolerance)


    totalrestrainedatoms=['N','CA','C']
    force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k", 1.0*kilocalories_per_mole/angstroms**2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    for i, atom_crd in enumerate(parm.positions):
        if parm.atoms[i].name in totalrestrainedatoms:
            force.addParticle(i, atom_crd.value_in_unit(u.nanometers))
    system.addForce(force)


    simulation = Simulation(topology, system, integrator, platform)#, platformProperties)
    simulation.context.setPositions(positions)
    if parm.box_vectors is not None:
        simulation.context.setPeriodicBoxVectors(*parm.box_vectors)


    simulation.reporters.append(dcdReporter)
    simulation.reporters.append(dataReporter)
    simulation.reporters.append(checkpointReporter)
    simulation.reporters.append(restartreporter)

    print("Starting Minimization .....")

    simulation.minimizeEnergy()
    simulation.saveState(min_state_path)
    simulation.saveCheckpoint(min_chk_path)

    simulation.context.setVelocitiesToTemperature(5*kelvin)
 
    T = 5

    print(f"Doing minimization and NVT for rep {rep}.....")
    for i in range(60):
      simulation.step(int(mdsteps/60))
      temperature = (T+(i*T))*kelvin
      #print("Temperature is: ", temperature)
      integrator.setTemperature(temperature)

    simulation.saveState(nvt_state_path)
    simulation.saveCheckpoint(nvt_chk_path)

