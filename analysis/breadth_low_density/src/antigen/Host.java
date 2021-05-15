package antigen;
/* A human individual that harbors viruses and immunity */

import java.util.*;
import java.io.*;

import org.tmatesoft.sqljet.core.SqlJetException;
import org.tmatesoft.sqljet.core.SqlJetTransactionMode;
import org.tmatesoft.sqljet.core.table.*;

public class Host {

	// fields
	private Virus infection;												
	private List<Phenotype> immuneHistory = new ArrayList<Phenotype>(0);
	private List<Double> infectionDates = new ArrayList<Double>(0);
	private List<Double> vaccinationDates = new ArrayList<Double>(0);
	private double birthDate;
	private IterableBitSet vaccinationHistory;
	private double lastVaccineDate;
	private double lastLastVaccineDate;
	
	// Naive host
	public Host(boolean vaccinationActive, double date) {
		if(vaccinationActive) {
			vaccinationHistory = new IterableBitSet(0);
		}
		birthDate = date;
	}
	
	// Immune host
	public Host(Phenotype immunity, boolean vaccinationActive) {
		assert immunity != null;
		addToImmuneHistory(immunity);
		if(vaccinationActive) {
			vaccinationHistory = new IterableBitSet(0);
		}
	}
	
	// Infected host
	public Host(Virus v, boolean vaccinationActive, double date) {
		infection = v;
		birthDate = date;
		if(vaccinationActive) {
			vaccinationHistory = new IterableBitSet(0);
		}
	}
	
//	public boolean equals(Virus v){
//		if(infection.equals(v)){
//			return true;
//		}else{
//			return false;
//		}
//	}
	
	public void addToImmuneHistory(Phenotype p) {
//		for(Phenotype pi : immuneHistory) {
//			if(pi.equals(p)) {
//				return;
//			}
//		}
		immuneHistory.add(p);
	}

	// infection methods
	public void reset(double date) {
		lastVaccineDate = 0;
		lastLastVaccineDate = 0;
		infection = null;
		infectionDates = new ArrayList<Double>(0);
		vaccinationDates = new ArrayList<Double>(0);
		immuneHistory = new ArrayList<Phenotype>(0);
		birthDate = date;
		if(vaccinationHistory != null) {
			vaccinationHistory = new IterableBitSet(0);
		}
	}
	
	public boolean isInfected() {
		boolean infected = false;
		if (infection != null) {
			infected = true;
		}
		return infected;
	}
	
	public Virus getInfection() {
		return infection;
	}
	
	public void infect(double time, Virus pV, int d) {
		Virus nV = new Virus(time, pV, d);
		infection = nV;
		infectionDates.add(time);
	}
	
	public void clearInfection() {
		Phenotype p = infection.getPhenotype();
		addToImmuneHistory(p);
		infection = null;
	}
	public int getHistoryLength() {
		return immuneHistory.size();
	}
	
	public void resetInfection(Phenotype p) {
		infection.setPhenotype(p);
		//System.err.println("Host received phenotype : "+infection.getPhenotype());
	}
	
	// make a new virus with the mutated phenotype
	public void mutate(boolean mut2D, double meanStep, double sdStep, double birth) {
		Virus mutV = infection.mutate(mut2D, meanStep, sdStep, birth);
		infection = mutV;
	}
	
	public void vaccinate(int vaccineId, double vaccineDate) {
		vaccinationHistory.set(vaccineId);
		lastLastVaccineDate = lastVaccineDate;
		lastVaccineDate = vaccineDate;
		vaccinationDates.add(vaccineDate);
	}
	
	public double getLastVaccineDate(){
		return lastVaccineDate;
	}
	
	public double getLastLastVaccineDate(){
		return lastLastVaccineDate;
	}
	
	// history methods
	
	public List<Phenotype> getHistory() {
		return immuneHistory;
	}	
	
	public IterableBitSet getVaccinationHistoryIdSet(Simulation sim) {
		return vaccinationHistory;
	}
	
	public List<Phenotype> getVaccinationHistory(Simulation sim) {
		if(vaccinationHistory == null) {
			return null;
		}
		
		List<Phenotype> vhList = new ArrayList<Phenotype>(vaccinationHistory.cardinality());
		for(int vaccineId : vaccinationHistory) {
			vhList.add(sim.getVaccine(vaccineId));
		}
		return vhList;
	}
	
	public void printImmuneHistory() {
		for (int i = 0; i < immuneHistory.size(); i++) {
			System.out.println(immuneHistory.get(i));
		}
	}

	public void printInfection(PrintStream stream) {
		if (infection != null) {
			stream.print(infection.getPhenotype());
		}
		else {
			stream.print("n");
		}
	}
	
	public void printHistory(PrintStream stream) {
		if (immuneHistory.size() > 0) {
			stream.print(immuneHistory.get(0));
			for (int i = 1; i < immuneHistory.size(); i++) {
				stream.print(";" + immuneHistory.get(i));
			}
		}
		else {
			stream.print("n");
		}
	}	
	
	public void printBirthDate(PrintStream stream) {
		stream.print(birthDate);
	}
	
	public void printInfectionDates(PrintStream stream) {
		if (infectionDates.size() > 0) {
			stream.print(infectionDates.get(0));
			for (int i = 1; i < infectionDates.size(); i++) {
				stream.print(";" + infectionDates.get(i));
			}
		}
		else {
			stream.print("n");
		}
	}
	
	public void printVaccinationDates(PrintStream stream) {
		if (vaccinationDates.size() > 0) {
			stream.print(vaccinationDates.get(0));
			for (int i = 1; i < vaccinationDates.size(); i++) {
				stream.print(";" + vaccinationDates.get(i));
			}
		}
		else {
			stream.print("n");
		}
	}	
	
	public void printInfection_newFormat(PrintStream stream) {
		if (infection != null) {
			stream.print(infection.getPhenotype());
		}
	}
	
	public void writeVirusesToSqlite(String demeName, int hostId, Simulation sim, SqlJetDb sampleDb){
		try {
			sampleDb.beginTransaction(SqlJetTransactionMode.WRITE);
			ISqlJetTable table = sampleDb.getTable("viruses");
			table.insert(hostId, sim.getDate(), infection.getPhenotype().getTraitA(), infection.getPhenotype().getTraitB());
			sampleDb.commit();
		}
		catch(SqlJetException e) {
			throw new RuntimeException(e);
		}
	}
	
	
	public void writeInfectedToSqlite(String demeName, int hostId, Simulation sim, SqlJetDb sampleDb){
		try {
			sampleDb.beginTransaction(SqlJetTransactionMode.WRITE);
			ISqlJetTable table = sampleDb.getTable("hosts");
			table.insert(demeName, hostId, sim.getDay(), infection.getPhenotype().getTraitA(), infection.getPhenotype().getTraitB(),"I");
			sampleDb.commit();
		}
		catch(SqlJetException e) {
			throw new RuntimeException(e);
		}
	}
	
	public void writeHostsToSqlite(String demeName, int hostId, Simulation sim, SqlJetDb sampleDb, double S) {
//		if(immuneHistory.size() > 0){
//			for (int i = 0; i < immuneHistory.size(); i++) {
//				try {
//					sampleDb.beginTransaction(SqlJetTransactionMode.WRITE);
//					ISqlJetTable table = sampleDb.getTable("hosts");
//					table.insert(demeName, hostId, sim.getDay(), immuneHistory.get(i).getTraitA(),immuneHistory.get(i).getTraitB(),"S");
//					sampleDb.commit();
//				}
//				catch(SqlJetException e) {
//					throw new RuntimeException(e);
//				}
//			}
//		}
		List<Phenotype> vH = getVaccinationHistory(sim);
		if(vH.size() > 0){
			for(int i = vH.size() - 1; i >= 0; i--){
				try {
					sampleDb.beginTransaction(SqlJetTransactionMode.WRITE);
					ISqlJetTable table = sampleDb.getTable("hosts");
					Phenotype v = vH.get(i);
					table.insert(demeName, hostId, sim.getDate(), v.getTraitA(), v.getTraitB(),S);
					sampleDb.commit();
				}
				catch(SqlJetException e) {
					throw new RuntimeException(e);
				}
			}
		}

//		if(immuneHistory.size() == 0 & vaccinationHistory.cardinality() == 0){
//			try {
//				sampleDb.beginTransaction(SqlJetTransactionMode.WRITE);
//				ISqlJetTable table = sampleDb.getTable("hosts");
//				table.insert(demeName, hostId, sim.getDate(), null, null,"none");
//				sampleDb.commit();
//			}
//			catch(SqlJetException e) {
//				throw new RuntimeException(e);
//			}
//		}
	}
		
	public String toString() {
		return Integer.toHexString(this.hashCode());
	}	
	
}