package antigen;
/* A population of host individuals */

import java.util.*;
import java.io.*;

import org.tmatesoft.sqljet.core.SqlJetException;
import org.tmatesoft.sqljet.core.SqlJetTransactionMode;
import org.tmatesoft.sqljet.core.table.*;


public class HostPopulation {
	private Simulation sim;
	private Parameters params;
	private Phenotype urImmunity;
	
	// fields
	private int deme;
	private String name;	
	private int cases;	
	public List<Host> susceptibles;
	public List<Host> infecteds;  
	private List<Host> recovereds;
	private List<Host> transcendentals;
	private double diversity;
	private double tmrca;
	private double netau;	
	private double serialInterval;
	private double antigenicDiversity;		
	
	private List<Phenotype> currentPhenotypes;
	
	private int newContacts;
	private int newRecoveries;

	private List<Integer> nVaccinatedS;
	private List<Integer> nVaccinatedI;
	
	public final int initialI;
	public final int initialR;
	public final int initialS;
	public final int initialT;
	
	private Double currentVaccinationRate;

	// construct population, using Virus v as initial infection
	public HostPopulation(Simulation sim, Parameters params, Phenotype urImmunity, Virus urVirus, int d) {
		this.sim = sim;
		this.params = params;
		this.urImmunity = urImmunity;
		
		// basic params
		deme = d;
		name = params.demeNames[deme];
		
		currentVaccinationRate = params.vaccinationRate[deme];
		
		// Calculate equilibrium
		assert params.swapDemography; // Code currently assumes this
		assert params.birthRate[d] == params.deathRate[d]; // Protect against human error specifying the wrong one
		
		double R0 = params.beta / (params.nu + params.deathRate[d]);
		assert R0 > 1.0;
		
		double initialFracS = 1.0 / R0;
		double initialFracI = params.deathRate[d] / params.beta * (R0 - 1.0);
		double riskAdjustment = Math.max(
			1.0 - params.smithConversion * Math.abs(params.initialTraitA),
			1.0 - params.homologousImmunity
		);
		assert riskAdjustment > 0.0;
		double initialFracR = (1.0 - initialFracS - initialFracI) / riskAdjustment;
		assert initialFracR > 0.0;
		
		if(params.startAtEquilibriumInfected[d]) {
			initialI = (int)Math.round(params.initialNs[d] * initialFracI);
			
		}
		else {
			initialI = params.initialIs[d];
		}
		System.err.printf("initial I: %d\n", initialI);
		
		if(params.startAtEquilibriumImmune[d]) {
			initialR = (int)Math.round(params.initialNs[d] * initialFracR);
		}
		else {
			initialR = (int)Math.round(params.initialPrR * params.initialNs[d]);
		}
		System.err.printf("initial R: %d\n", initialR);
		
		assert initialI + initialR <= params.initialNs[d];
		initialS = params.initialNs[d] - initialI - initialR; 
		System.err.printf("initial S: %d\n", initialS);
		
		if (params.transcendental) {
			initialT = (int) ((double) params.initialNs[deme] * params.initialPrR);
		}
		else {
			initialT = 0;
		}
		System.err.printf("initial T: %d\n", initialT);
		
		// Initialize susceptible individuals
		susceptibles = new ArrayList<Host>(params.initialNs[d]);
		for(int i = 0; i < initialS; i++) {
			Host h = new Host(params.vaccinate, sim.getDate());
			susceptibles.add(h);
		}
		
		// Initialize immune (or transcendental) individuals
		transcendentals = new ArrayList<Host>(initialT * 10);
		for(int i = 0; i < initialR; i++) {
			Host h = new Host(urImmunity, params.vaccinate);
			if (params.transcendental) {
				transcendentals.add(h);
			}
			else {
				susceptibles.add(h);
			}
		}

		// Initialize infected individuals
		infecteds = new ArrayList<Host>(initialI * 10);
		for (int i = 0; i < initialI; i++) {
			Virus v = new Virus(0.0, urVirus, deme);
			Host h = new Host(v, params.vaccinate, sim.getDate());
			infecteds.add(h);
		}
		
		currentPhenotypes = new ArrayList<Phenotype>();
		
		nVaccinatedS = new ArrayList<Integer>();
		nVaccinatedI = new ArrayList<Integer>();
	}
	
	// accessors
	public int getUnexposed(){
		int unexposed = 0;
		//TODO: reimplement if needed
//		for(int i =0; i<susceptibles.size(); i++){
//			Host h = susceptibles.get(i);
//			if(h.getHistory().size() == 0){
//				unexposed++;
//			}
//		}
		return unexposed;
	}
	
	public List<Phenotype> getCurrentPhenotypes(){
		return currentPhenotypes;
	}
	
	public int getN() {
		return susceptibles.size() + infecteds.size() + transcendentals.size();
	}
	public int getS() {
		return susceptibles.size();
	}
	public int getI() {
		return infecteds.size();
	}
	public int getT() {
		return transcendentals.size();
	}
	public double getPrS() {
		return (double) getS() / (double) getN();
	}
	public double getPrI() {
		return (double) getI() / (double) getN();
	}
	public double getPrT() {
		return (double) getT() / (double) getN();
	}
	public int getRandomN() {
		return Random.nextInt(0,getN()-1);
	}
	public int getRandomS() {
		return Random.nextInt(0,getS()-1);
	}
	public int getRandomI() {
		return Random.nextInt(0,getI()-1);
	}
	public int getRandomT() {
		return Random.nextInt(0,getT()-1);
	}
	
	public Host getRandomHost() {
		int i = Random.nextInt(0, getN() - 1);
		if(i <= getS()) {
			return getRandomHostS();
		}
		else if(i > getS() && i <= getS() + getI()){
			return getRandomHostI();
		}
		return getRandomHostT();
	}
	
	public Host getRandomHostS() {
		int index = Random.nextInt(0,getS()-1);
		return susceptibles.get(index);
	}
	public Host getRandomHostI() {
		Host h = null;
		if (getI() > 0) {
			int index = Random.nextInt(0,getI()-1);
			h = infecteds.get(index);
		}
		return h;
	}
	public Host getRandomHostT() {
		Host h = null;
		if (getT() > 0) {
			int index = Random.nextInt(0,getT()-1);
			h = infecteds.get(index);
		}
		return h;
	}
	
	public Virus getRandomInfection() {
		Virus v = null;
		Host h = getRandomHostI();
		if (h != null) {
			v = h.getInfection();
		}
		return v;
	}	
	
	public void resetCases() {
		cases = 0;
	}
	public int getCases() {
		return cases;
	}	

	public double getDiversity() {
		return diversity;
	}		
	
	public double getNetau() {
		return netau;
	}	
	
	public double getTmrca() {
		return tmrca;
	}	
	
	public double getSerialInterval() {
		return serialInterval;	
	}		
	
	public double getAntigenicDiversity() {
		return antigenicDiversity;
	}			
	
	public double getVaccinationRate() {
		return currentVaccinationRate;
	}
	
	public void removeSusceptible(int i) {
		int lastIndex = getS() - 1;
		Host lastHost = susceptibles.get(lastIndex);
		susceptibles.set(i,lastHost);
		susceptibles.remove(lastIndex);
	}	
	public void removeInfected(int i) {
		int lastIndex = getI() - 1;
		Host lastHost = infecteds.get(lastIndex);
		infecteds.set(i,lastHost);
		infecteds.remove(lastIndex);
	}
	public void removeTranscendental(int i) {
		int lastIndex = getT() - 1;
		Host lastHost = transcendentals.get(lastIndex);
		transcendentals.set(i,lastHost);
		transcendentals.remove(lastIndex);	
	}
	
	public void stepForward(double day) {
//		if(currentPhenotypes.size()>1){
//			System.err.println(name+" " + sim.day+":"+currentPhenotypes.get(1));
//		}
		
	//	resetCases();
		if (params.swapDemography) {
			swap();
		} else {
			assert false;
//			grow();
//			decline();
		}
		
		//updateAges();
		recordContacts();
		recordRecoveries();
		distributeContacts();
		distributeRecoveries();	
		loseImmunity();
		
		mutate();
		sample();

	
	}
	
	// draw a Poisson distributed number of births and add these hosts to the end of the population list
	// TODO: reactivate this if necessary
//	public void grow() {
//		double totalBirthRate = getN() * params.birthRate[deme] * params.deltaT;
//		int births = Random.nextPoisson(totalBirthRate);
//		for (int i = 0; i < births; i++) {
//			Host h = new Host(params.initialPrR, urImmunity, params.vaccinate);
//			susceptibles.add(h);
//		}
//	}
//	
//	// draw a Poisson distributed number of deaths and remove random hosts from the population list
//	public void decline() {
//		// deaths in susceptible class
//		double totalDeathRate = getS() * params.deathRate[deme] * params.deltaT;
//		int deaths = Random.nextPoisson(totalDeathRate);
//		for (int i = 0; i < deaths; i++) {
//			if (getS()>0) {
//				int index = getRandomS();
//				removeSusceptible(index);
//			}
//		}		
//		// deaths in infectious class		
//		totalDeathRate = getI() * params.deathRate[deme] * params.deltaT;
//		deaths = Random.nextPoisson(totalDeathRate);
//		for (int i = 0; i < deaths; i++) {
//			if (getI()>0) {
//				int index = getRandomI();
//				removeInfected(index);
//			}
//		}
// 		// deaths in recovered class
// 		totalDeathRate = getR() * params.deathRate[deme] * params.deltaT;
// 		deaths = Random.nextPoisson(totalDeathRate);
// 		for (int i = 0; i < deaths; i++) {
// 				if (getR()>0) {
// 				int index = getRandomR();
// 				removeRecovered(index);
// 			}
// 		}
//	}
	
	// draw a Poisson distributed number of births and reset these individuals
	public void swap() {
		// draw random individuals from susceptible class
		double totalBirthRate = getS() * params.birthRate[deme] * params.deltaT;
		int births = Random.nextPoisson(totalBirthRate);
		for (int i = 0; i < births; i++) {
			if (getS()>0) {
				int index = getRandomS();
				Host h = susceptibles.get(index);
				if(params.vaccinate) {
					for(int vaccineId : h.getVaccinationHistoryIdSet(sim)) {
						nVaccinatedS.set(vaccineId, nVaccinatedS.get(vaccineId) - 1);
					}
				}
				
				h.reset(sim.getDate());
			}
		}		
		// draw random individuals from infected class
		totalBirthRate = getI() * params.birthRate[deme] * params.deltaT;
		births = Random.nextPoisson(totalBirthRate);
		for (int i = 0; i < births; i++) {
			if (getI()>0) {
				int index = getRandomI();
				Host h = infecteds.get(index);
				if(params.vaccinate) {
					for(int vaccineId : h.getVaccinationHistoryIdSet(sim)) {
						nVaccinatedI.set(vaccineId, nVaccinatedI.get(vaccineId) - 1);
					}
				}
				h.reset(sim.getDate());
				removeInfected(index);
				susceptibles.add(h);
			}
		}
		// draw random individuals from recovered class
		totalBirthRate = getT() * params.birthRate[deme] * params.deltaT;
		births = Random.nextPoisson(totalBirthRate);
		for (int i = 0; i < births; i++) {
			if (getT()>0) {
				int index = getRandomT();
				Host h = transcendentals.get(index);
				h.reset(sim.getDate());
				removeTranscendental(index);
				susceptibles.add(h);
			}
		}
		
	}
	
	Set<Integer> chooseRandomIndices(int size, int count) {
		assert(count <= size);
		Set<Integer> indices = new HashSet<>();
		for(int i = 0; i < count; i++) {
			int index;
			do {
				index = Random.nextInt(0, size - 1);
			} while(indices.contains(index));
			indices.add(index);
		}
		return indices;
	}
	
	public double getFrequency(Phenotype p, List<Phenotype> pList){
		double pCounts = (double) Collections.frequency(pList, p);
		return (pCounts / pList.size());
//		double freq = (double) Collections.frequency(infecteds, v);
//		return freq / infecteds.size();
	}
	
	public void updateStrains(){
		if(infecteds.size() != 0 ){
			currentPhenotypes = new ArrayList<Phenotype>();
			Set<Integer> sampleIndices = chooseRandomIndices(infecteds.size(), (int) Math.round(infecteds.size()*params.strainSampleRate));
			for(int index : sampleIndices) {
				currentPhenotypes.add(infecteds.get(index).getInfection().getPhenotype());
			}
		}
//		System.err.println("nviruses = "+ currentPhenotypes.size());
//		Phenotype p = currentPhenotypes.get(0);
//		System.err.println(sim.getDay()+" "+p+" is at: "+getFrequency(p, currentPhenotypes));
	}
	
	public boolean resetStrains(SqlJetDb readDb){
		if(sim.getDay() == 0){
			return(true);
		} else {
			int nInfections = getI();
			if (nInfections > 0) {
				int index = 0;
				try {
					readDb.beginTransaction(SqlJetTransactionMode.READ_ONLY);
					ISqlJetTable phenotypes = readDb.getTable("phenotypes");
					ISqlJetCursor c = phenotypes.lookup("day", sim.getDay());
					if (c.getRowCount() > 0) {
						do {
							Phenotype p = new Phenotype(
									(double) c.getFloat("ag1"),
									(double) c.getFloat("ag2"));
							int niStrain = (int) Math.round(c.getFloat("freq")
									* nInfections);

							for (int j = 0; j < niStrain; j++) {
								if (index < nInfections) {
									infecteds.get(index).resetInfection(p);
									index++;
								}
							}
						} while (c.next());
						readDb.commit();
						return (true);
					} else {
						readDb.commit();
						return (false);
					}
				} catch (SqlJetException e) {
					throw new RuntimeException(e);
				}
			}
			return (false);
		}
	}
	
	public void resetStrains(List<Phenotype> reference){
		int nInfections = getI();
		
		if(nInfections > 0){
			int index = 0;
			for (Phenotype iStrain : reference) {
				//System.err.println(iStrain);
				if(iStrain.getDate() == sim.getDate()){
					// System.err.println("resetting:"+iStrain);
					// System.err.println(sim.day +" "+ iStrain+" freq: "+ Math.round(getFrequency(iStrain, reference)*nInfections)+" "+getI());
					int niStrain = (int) Math.round(getFrequency(iStrain, reference) * nInfections);
					
					//System.err.println(sim.day+":"+i+":"+niStrain+":"+nInfections);
					for(int j = 0; j < niStrain; j++){
						if(index < nInfections){
							
							infecteds.get(index).resetInfection(iStrain);
							index ++;
						}
					}
				}
			}
		}
	}
	
	public void updateVaccinationRate(){
		currentVaccinationRate = 0.0;
		//currentVaccinationRate = params.allVaccinationRates[Random.nextInt(0,params.allVaccinationRates.length-1)] * params.deltaT * 365.0/ params.vaccineWindow;
		System.err.println("Updating rate :" + currentVaccinationRate);
	}
	
	private void vaccinate(int vaccineId, List<Host> hosts, List<Integer> counts, double vaccinationRate) {
		if(hosts.size() == 0) {
			return;
		}
		
//		System.err.printf("Vaccination rate: %f\n", vaccinationRate);
		int vacCount = Math.min(hosts.size(), Random.nextPoisson(hosts.size() * vaccinationRate));
//		System.err.printf("Vaccination count: %d\n", vacCount);
		Set<Integer> vacIndices = chooseRandomIndices(hosts.size(), vacCount); 
		for(int index : vacIndices) {
			hosts.get(index).vaccinate(vaccineId, sim.getDate());
		}
		if(vaccineId == counts.size()) {
			counts.add(vacCount);	
		}
		else if(vaccineId == counts.size() - 1) {
			counts.set(vaccineId, counts.get(vaccineId) + vacCount);
		}
		else {
			assert false;
		}
	}
	
	// vaccinate a Poisson-distributed number of hosts (in all states)
	public void vaccinate(int vaccineId) {
		double vaccinationRate = currentVaccinationRate * params.deltaT * 365.0 / params.vaccineWindow;
		if(params.varyVaccinationRate){
			vaccinationRate = currentVaccinationRate * params.deltaT * 365.0/ params.vaccineWindow;
		}
		
		vaccinate(vaccineId, susceptibles, nVaccinatedS, vaccinationRate);
		vaccinate(vaccineId, infecteds, nVaccinatedI, vaccinationRate);
	}
	
	// draw a Poisson distributed number of contacts
	public void recordContacts() {
		// each infected makes I->S contacts on a per-day rate of beta * S/N
		double totalContactRate = getI() * getPrS() * params.beta * getSeasonality() * params.deltaT;
		newContacts = Random.nextPoisson(totalContactRate);			
	}

	// move from S->I following number of new contacts
	public void distributeContacts() {
		
		for (int i = 0; i < newContacts; i++) {
			if (getS()>0 && getI()>0) {
		
				// get indices and objects
				int index = getRandomI();
				int sndex = getRandomS();			
				Host iH = infecteds.get(index);	
				Host sH = susceptibles.get(sndex);	
//				if(params.useReferenceStrains){
//					index = sim.getDeme(0).getRandomI();
//					iH = sim.getDeme(0).infecteds.get(index);
//				}
				Virus v = iH.getInfection();
								
				// attempt infection
				Phenotype p = v.getPhenotype();		
				List<Phenotype> immuneHistory = sH.getHistory(); 
				List<Phenotype> vaccinationHistory = sH.getVaccinationHistory(sim);
				double lastVaccineDate = sH.getLastVaccineDate();
				double date = sim.getDate();
				double chanceOfSuccess = p.riskOfInfection(immuneHistory, vaccinationHistory, lastVaccineDate, date,
						params.smithConversion, params.homologousImmunity, params.vaccineImmuneBreadth
				);
				if (Random.nextBoolean(chanceOfSuccess)) {
					sH.infect(sim.getDate(), v,deme);
					
					if(params.vaccinate) {
						for(int vaccineId : sH.getVaccinationHistoryIdSet(sim)) {
							nVaccinatedS.set(vaccineId, nVaccinatedS.get(vaccineId) - 1);
							nVaccinatedI.set(vaccineId, nVaccinatedI.get(vaccineId) + 1);
						}
					}
					
					removeSusceptible(sndex);
					infecteds.add(sH);
					cases++;
				}
			
			}
		}		
		
	}
	
	// draw a Poisson distributed number of contacts and move from S->I based upon this
	// this deme is susceptibles and other deme is infecteds
	public void betweenDemeContact(HostPopulation hp, double day) {

		// each infected makes I->S contacts on a per-day rate of beta * S/N
		// double totalContactRate = hp.getI() * getPrS() * params.beta * params.betweenDemePro * params.demeBaselines[deme] * params.deltaT;//makes migration aseasonal
		// System.err.println("Contacting from: "+ hp.deme + " to " + this.deme + " at " + params.contactMatrix[hp.deme][this.deme]);
		
		double totalContactRate = hp.getI() * getPrS() * params.beta * params.betweenDemePro * params.contactMatrix[hp.deme][this.deme] *  getSeasonality() * params.deltaT;

		int contacts = Random.nextPoisson(totalContactRate);
		for (int i = 0; i < contacts; i++) {
			if (getS()>0 && hp.getI()>0) {
		
				// get indices and objects
				Host iH = hp.getRandomHostI();
				int sndex = getRandomS();
				Host sH = susceptibles.get(sndex);			
				Virus v = iH.getInfection();
				
				// attempt infection
				Phenotype p = v.getPhenotype();
				List<Phenotype> immuneHistory = sH.getHistory();
				List<Phenotype> vaccinationHistory = sH.getVaccinationHistory(sim);
				double lastVaccineDate = sH.getLastVaccineDate();
				double date = sim.getDate();
				
				double chanceOfSuccess = p.riskOfInfection(
					immuneHistory, vaccinationHistory, lastVaccineDate, date,
					params.smithConversion,
					params.homologousImmunity,
					params.vaccineImmuneBreadth
				);
				if (Random.nextBoolean(chanceOfSuccess)) {
					sH.infect(sim.getDate(), v,deme);
					if(params.vaccinate){
						for(int vaccineId : sH.getVaccinationHistoryIdSet(sim)) {
							nVaccinatedS.set(vaccineId, nVaccinatedS.get(vaccineId) - 1);
							nVaccinatedI.set(vaccineId, nVaccinatedI.get(vaccineId) + 1);
						}
					}
					removeSusceptible(sndex);
					infecteds.add(sH);
					cases++;
				}
			
			}
		}		
		
	}	
	
	// draw a Poisson distributed number of recoveries
	public void recordRecoveries() {	
		// each infected recovers at a per-day rate of nu
		double totalRecoveryRate = getI() * params.nu * params.deltaT;
		newRecoveries = Random.nextPoisson(totalRecoveryRate);	
	}
	
	// move from I->S following number of recoveries
	public void distributeRecoveries() {

		for (int i = 0; i < newRecoveries; i++) {
			if (getI()>0) {
				int index = getRandomI();
				Host h = infecteds.get(index);
				h.clearInfection();
				removeInfected(index);
				
				if (params.transcendental) {
					transcendentals.add(h);
				} else {
					susceptibles.add(h);
				}
				
				if(params.vaccinate) {
					for(int vaccineId : h.getVaccinationHistoryIdSet(sim)) {
						nVaccinatedI.set(vaccineId, nVaccinatedI.get(vaccineId) - 1);
						nVaccinatedS.set(vaccineId, nVaccinatedS.get(vaccineId) + 1);
					}
				}
			}
		}			
	}
	
	// draw a Poisson distributed number of T->S 
	public void loseImmunity() {
		// each recovered regains immunity at a per-day rate
		double totalReturnRate = getT() * params.immunityLoss * params.deltaT;
		int returns = Random.nextPoisson(totalReturnRate);
		for (int i = 0; i < returns; i++) {
			if (getT()>0) {
				int index = getRandomT();
				Host h = transcendentals.get(index);
				removeTranscendental(index);
				susceptibles.add(h);
			}
		}			
	}	
	
	
	// draw a Poisson distributed number of mutations and mutate based upon this
	// mutate should not impact other Virus's Phenotypes through reference
	public void mutate() {
		// each infected mutates at a per-day rate of mu
		double totalMutationRate = getI() * params.muPhenotype[deme] * params.deltaT;
		int mutations = Random.nextPoisson(totalMutationRate);
		for (int i = 0; i < mutations; i++) {
			if (getI()>0) {
				int index = getRandomI();
				Host h = infecteds.get(index);
				h.mutate(params.mut2D, params.meanStep, params.sdStep, sim.getDate());
			}
		}			
	}	
	
	// draw a Poisson distributed number of samples and add them to the VirusSample
	// only sample after burnin is completed
	public void sample() {
		if (getI()>0 && sim.getDay() >= params.burnin) {
		
			double totalSamplingRate = params.tipSamplingRate * params.deltaT;
			if (params.tipSamplingProportional) {
				totalSamplingRate *= getI();
			} 
			
			int samples = Random.nextPoisson(totalSamplingRate);
			for (int i = 0; i < samples; i++) {
				int index = getRandomI();
				Host h = infecteds.get(index);
				Virus v = h.getInfection();
				VirusTree.add(v);
			}
		}
	}
		
	// through current infected population assigning ancestry as trunk
	public void makeTrunk() {
		for (int i = 0; i < getI(); i++) {
			Host h = infecteds.get(i);
			Virus v = h.getInfection();
			v.makeTrunk();
			while (v.getParent() != null) {
				v = v.getParent();
				if (v.isTrunk()) {
					break;
				} else {
					v.makeTrunk();
				}
			}
		}
	}	
	
	public void updateDiversity() {

		diversity = 0.0;
		tmrca = 0.0;
		antigenicDiversity = 0.0;		
		netau = 0.0;
		serialInterval = 0.0;
		
		if (getI()>1) { 
		
			double coalCount = 0.0;	
			double coalOpp = 0.0;
			double coalWindow = params.netauWindow / 365.0;
			int sampleCount = params.diversitySamplingCount;
			
			for (int i = 0; i < sampleCount; i++) {
				Virus vA = getRandomInfection();
				Virus vB = getRandomInfection();
				if (vA != null && vB != null) {
					double dist = vA.distance(vB);
					diversity += dist;
					if (dist > tmrca) {
						tmrca = dist;
					}
					antigenicDiversity += vA.antigenicDistance(vB);
					coalOpp += coalWindow;
					coalCount += vA.coalescence(vB, coalWindow);
					serialInterval += vA.serialInterval();
				}
			}	
		
			diversity /= (double) sampleCount;
			tmrca /= 2.0;
			antigenicDiversity /= (double) sampleCount;		
			netau = coalOpp / coalCount;
			serialInterval /= (double) sampleCount;
		
		}
		
	}	
		
	public void printState(PrintStream stream) {
		updateDiversity();
		stream.printf("\t%.4f\t%.4f\t%.4f\t%.5f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d", getDiversity(), getTmrca(), getNetau(), getSerialInterval(), getAntigenicDiversity(), getN(), getS(), getI(), getT(), getCases(), getUnexposed());
	}	
	
	public void printHeader(PrintStream stream) {
		stream.printf("\t%sDiversity\t%sTmrca\t%sNetau\t%sSerialInterval\t%sAntigenicDiversity\t%sN\t%sS\t%sI\t%sT\t%sCases\t%sUnexposed", name, name, name, name, name, name, name, name, name, name, name);
	}
	
	
	public void printViruses(PrintStream stream){
		stream.print("deme,ag1,ag2\n");
		// step through infecteds and print
		for (int i = 0; i < getI(); i++) {
			Host h = infecteds.get(i);
			stream.print(deme + ",");
			h.printInfection(stream);
			stream.print("\n");
		}
	}
	
	public void writeHostsToSqlite(SqlJetDb sampleDb){

		// step through susceptibles and print
		for (int i = 0; i < params.writeSampleRateS*getS(); i++) {
			Host h = getRandomHostS();
			h.writeHostsToSqlite(name, i, sim, sampleDb, getS());
		}
		
		// step through infecteds and print
//		if(infecteds.size() != 0 ){
//			for (int i = 0; i < params.writeSampleRateI*getI(); i++) {
//				Host h = infecteds.get(i);
//				h.writeInfectedToSqlite(name, i, sim, sampleDb);
//			}	
//		}
		
	}
	
	public void writeVirusesToSqlite(SqlJetDb sampleDb){

		// step through infecteds and print
		if(infecteds.size() != 0 ){
			for (int i = 0; i < params.writeSampleRateI*getI(); i++) {
				Host h = infecteds.get(i);
				h.writeVirusesToSqlite(name, i, sim, sampleDb);
			}	
		}
		
	}
	
	public void writePhenotypesToSqlite(SqlJetDb sampleDb, PrintStream stream){
		Set<Phenotype> currentPSet = new HashSet<Phenotype>(currentPhenotypes);
		try{
			sampleDb.beginTransaction(SqlJetTransactionMode.WRITE);
			for (Phenotype iStrain : currentPSet) {
				
				// System.err.println(sim.day +" "+ iStrain+" freq: "+ Math.round(getFrequency(iStrain, reference)*nInfections)+" "+getI());
				double fiStrain = getFrequency(iStrain, currentPhenotypes);

				ISqlJetTable table = sampleDb.getTable("phenotypes");
				table.insert(sim.getDay(), iStrain.getTraitA(), iStrain.getTraitB(), fiStrain);
			}
			sampleDb.commit();
		}
		catch(SqlJetException e) {
			throw new RuntimeException(e);
		}

	}
	
	public void printLongitudinal(PrintStream stream){
		double sampleRateS = 0.05;
		double sampleRateI = 0.05;
		Set<Integer> sampleIndices = chooseRandomIndices(susceptibles.size(), (int) Math.round(susceptibles.size()*sampleRateS));

		// step through susceptibles and print
		for (int index : sampleIndices) {
			Host h = susceptibles.get(index);
			stream.print(index + ":");
			h.printInfectionDates(stream);
			stream.print(":");
			h.printVaccinationDates(stream);
			stream.print(":");
			h.printBirthDate(stream);
			stream.print("\n");
		}
		
		sampleIndices = chooseRandomIndices(infecteds.size(), (int) Math.round(infecteds.size()*sampleRateI));
		// step through infecteds and print
		for (int index : sampleIndices) {
			Host h = infecteds.get(index);
			stream.print((index+getS()) + ":");
			h.printInfectionDates(stream);
			stream.print(":");
			h.printVaccinationDates(stream);
			stream.print(":");
			h.printBirthDate(stream);
			stream.print("\n");
		}
	}
	
	public void printSomeHosts(PrintStream stream){
		Set<Integer> sampleIndices = chooseRandomIndices(susceptibles.size(), (int) Math.round(susceptibles.size()*params.writeSampleRateS));
		
		for (int index : sampleIndices) {
			Host h = susceptibles.get(index);
			stream.printf("%.4f\t%d\t%d\t%f\t%f\t%f\t%d\t%d\t%f\t%f\t%f\t%f\n",
					sim.getDate(), 
					(int) Math.min(h.getVaccinationHistory(sim).size(),1), 
					0,
					h.getLastVaccineDate(),
					h.getLastLastVaccineDate(),
					currentVaccinationRate, 
					(int) h.getVaccinationHistory(sim).size(), 
					(int) h.getHistoryLength(),
					-100.0,
					-100.0,
					sim.getLastVaccineStrain().getTraitA(),
					sim.getLastVaccineStrain().getTraitB());
		}
		
		sampleIndices = chooseRandomIndices(infecteds.size(), (int) Math.round(infecteds.size()*params.writeSampleRateI));
		for (int index : sampleIndices) {
			Host h = infecteds.get(index);
			stream.printf("%.4f\t%d\t%d\t%f\t%f\t%f\t%d\t%d\t%f\t%f\t%f\t%f\n",
					sim.getDate(), 
					(int) Math.min(h.getVaccinationHistory(sim).size(),1), 
					1,
					h.getLastVaccineDate(),
					h.getLastLastVaccineDate(),
					currentVaccinationRate,
					(int) h.getVaccinationHistory(sim).size(), 
					(int) h.getHistoryLength(),
					h.getInfection().getPhenotype().getTraitA(),
					h.getInfection().getPhenotype().getTraitB(),
					sim.getLastVaccineStrain().getTraitA(),
					sim.getLastVaccineStrain().getTraitB());
		}
		
	}

	
	public void printHostPopulation(PrintStream stream) {
		
		// step through susceptibles and print
		stream.print("Susceptibles\n");
		for (int i = 0; i < getS(); i++) {
			Host h = getRandomHostS();
			stream.print(deme + ":");
			h.printInfection(stream);
			stream.print(":");
			h.printHistory(stream);
			stream.print(":");
		}
		
		// step through infecteds and print
		stream.print("Infecteds\n");
		for (int i = 0; i < getI(); i++) {
			Host h = infecteds.get(i);
			stream.print(deme + ":");
			h.printInfection(stream);
			stream.print(":");
			h.printHistory(stream);
			stream.print(":");
		}
		
		// step through recovereds and print
		for (int i = 0; i < getT(); i++) {
			Host h = recovereds.get(i);
			stream.print(deme + ":");
			h.printInfection(stream);
			stream.print(":");
			h.printHistory(stream);
			stream.println();
		}		
		
	}
	
	public Double getSeasonality() {
		double baseline = params.demeBaselines[deme];
		double amplitude = params.demeAmplitudes[deme];
		double offset = params.demeOffsets[deme];
		double beta = baseline + amplitude * Math.cos(2*Math.PI*sim.getDate() + 2*Math.PI*offset);
		return beta;
	}
	
	public int getNVaccinatedS(int vaccineId) {
		return nVaccinatedS.get(vaccineId);
	}
	
	public int getNVaccinatedI(int vaccineId) {
		return nVaccinatedI.get(vaccineId);
	}
				
}