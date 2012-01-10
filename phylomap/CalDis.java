package phylomap;
import java.util.*;
import java.io.*;


public class CalDis {

	/**
	 * @param args
	 */
	
public int seqdiff_old(String s1,String s2){
		
		int score=0;
		int kscore=0;
		int lkscore=0;
		int adscore=0;
		if(s1.length()!=s2.length()){
			return -1;
		}else{
			
			for(int i=0;i<s1.length();i++){
				String c1=s1.substring(i,i+1);
				String c2=s2.substring(i,i+1);
				if(c1.equals(c2)){
					
				}else{
					if(c1.equals("-")||c2.equals("-")){
						if((i>=40&&i<=200) ){
							score++;
						}else if(i<40){
							lkscore++;
						}
						else if(i<230){
							kscore++;
						}else if(i>=230){
							adscore++;
						}
					}else{
						score++;
					}
				}
				
			}
			if(kscore==0){
				return score+(lkscore/8)+adscore;
			}else{
				return score+((kscore+lkscore+adscore)/8);
			}
		}
		
		
	}
	

	// compare 2 id length seq to find how many diff they have
	public int seqdiff(String s1,String s2){
		s1=s1.toUpperCase();
		s2=s2.toUpperCase();
		int score=0;
		if(s1.length()!=s2.length()){
			System.out.println("Fatal error, seq not alinged!");
			System.out.println(s1.length()+":"+s2.length());
			return -1;
		}else{
			for(int i=0;i<s1.length();i++){
				String c1=s1.substring(i,i+1);
				String c2=s2.substring(i,i+1);
				if(!c1.equals(c2)){
					score++;
				}
			}
			return score;
		}
	}
	
	public int seqdiffEFT(String s1,String s2){
		
		int score=0;
		int kscore=0;

		if(s1.length()!=s2.length()){
			return -1;
		}else{
			
			for(int i=0;i<s1.length();i++){
				String c1=s1.substring(i,i+1);
				String c2=s2.substring(i,i+1);
				if(c1.equals(c2)){
					
				}else{
					if(c1.equals("-")||c2.equals("-")){
						if(i<=120){
							score++;
						}
						else{
							kscore++;
						}
					}else{
						score++;
					}
				}
				
			}
			return score+(kscore/5);
		}
		
		
	}

	public void calDisofSeqList(List l, File fout, String seqx, String seqy) throws Exception{
		StringBuffer sb=new StringBuffer("");
		for(int i=0;i<l.size();i++){
			String seq=(String)l.get(i);
			int x=this.seqdiff(seqx, seq);
			int y=this.seqdiff(seqy, seq);
			double xd=x+Math.random();
			if( xd<=0){
				xd=0.1;
			}
			double yd=y+Math.random();
			if(yd<=0){
				yd=0.1;
			}
			if(xd<100 && yd<100){
				sb.append(""+xd+"	"+yd+"\n");
			}
			System.out.println(i);
		}
		FileWriter fw4=new FileWriter(fout);
		fw4.write(sb.toString());
		fw4.close();
	}
	
	public void discretionMap(List l1, List l2, File fout) throws Exception{
		FileWriter fw=new FileWriter(fout);
		for(int i=0;i<l1.size();i++){
			String s=(String)l1.get(i);
			//inner 
			double x=this.avg_point_to_group_distance(s, l1);
			x+=Math.random();
			//outer 
			double y=this.avg_point_to_group_distance(s, l2);
			y+=Math.random();
			
			fw.write(""+x+"	"+y+"\n");
			System.out.println("list 1:"+i);
		}
		
		for(int i=0;i<l2.size();i++){
			String s=(String)l2.get(i);
			//outer
			double x=this.avg_point_to_group_distance(s, l1);
			x+=Math.random();
			//inner 
			double y=this.avg_point_to_group_distance(s, l2);
			y+=Math.random();
			
			fw.write(""+x+"	"+y+"\n");
			System.out.println("list 2:"+i);
		}
		
		fw.close();
	}
	
	public double avg_point_to_group_distance(String seq_point, List group){
		
		int sum=0;
		for(int i=0;i<group.size();i++){
			String s=(String)group.get(i);
			sum+=this.seqdiff(s, seq_point);
		}
		return ((double)sum)/((double)group.size());
	}
	
	public HashMap genSmithDisMatrix(){
		HashMap hm=new HashMap();
		//assign dis 1 
		hm.put("DE", new Integer(1));
		hm.put("ED", new Integer(1));
		hm.put("KH", new Integer(1));
		hm.put("HK", new Integer(1));
		hm.put("KR", new Integer(1));
		hm.put("RK", new Integer(1));
		hm.put("RH", new Integer(1));
		hm.put("HR", new Integer(1));
		hm.put("NQ", new Integer(1));
		hm.put("QN", new Integer(1));
		hm.put("ST", new Integer(1));
		hm.put("TS", new Integer(1));
		hm.put("AG", new Integer(1));
		hm.put("GA", new Integer(1));
		hm.put("IL", new Integer(1));
		hm.put("LI", new Integer(1));
		hm.put("IV", new Integer(1));
		hm.put("VI", new Integer(1));
		hm.put("LV", new Integer(1));
		hm.put("VL", new Integer(1));
		hm.put("FW", new Integer(1));
		hm.put("WF", new Integer(1));
		hm.put("FY", new Integer(1));
		hm.put("YF", new Integer(1));
		hm.put("WY", new Integer(1));
		hm.put("YW", new Integer(1));
		hm=this.fillHM(hm, "DEKRHNQST", new Integer(2));
		hm=this.fillHM(hm, "ILVFWYCM", new Integer(2));
		hm=this.fillHM(hm, "DEKRHNQSTILVFWYCMAGPX-", new Integer(3));
		return hm;
	}
	
	public HashMap fillHM(HashMap hm, String s, Integer score){
		
		for(int i=0;i<s.length();i++){
			String s1=s.substring(i,i+1);
			for(int j=i+1;j<s.length();j++){
				String s2=s.substring(j,j+1);
				//System.out.println(s1+s2);
				Object o1=hm.get(s1+s2);
				Object o2=hm.get(s2+s1);
				if(o1==null && o2==null){
					hm.put(s1+s2, score);
					hm.put(s2+s1, score);
				}
			}
			
			
		}
		
		return hm;
	}
	
	public double seqdiffSmith(String s1,String s2,HashMap hm){
		//HashMap hm=this.genSmithDisMatrix();
		s1=s1.toUpperCase();
		s2=s2.toUpperCase();
		double score=0;
		if(s1.length()!=s2.length()){
			System.out.println("Fatal error, seq not alinged!");
			System.out.println(s1.length()+":"+s2.length());
			System.out.println(s1);
			System.out.println(s2);
			return -1;
		}else{
			for(int i=0;i<s1.length();i++){
				String c1=s1.substring(i,i+1);
				String c2=s2.substring(i,i+1);
				if(!c1.equals(c2)){
					Integer Iscore=(Integer)hm.get(c1+c2);
					if(Iscore!=null){
						score+=Iscore.intValue();
					}else{
						score+=3;
					}
				}
			}
			//score=Math.sqrt(score/((double)s1.length()));
			score=Math.sqrt(score);
			return score;
		}
	}
	
	public double p_distance(String s1,String s2){
		s1=s1.toUpperCase();
		s2=s2.toUpperCase();
		int score=0;
		int cnt=0;
		if(s1.length()!=s2.length()){
			System.out.println("Fatal error, seq not alinged!");
			System.out.println(s1.length()+":"+s2.length());
			return -1;
		}else{
			for(int i=0;i<s1.length();i++){
				String c1=s1.substring(i,i+1);
				String c2=s2.substring(i,i+1);
				if("-".equals(c1) || "-".equals(c2)){
					cnt++;
				}else {
					cnt++;
					if(!c1.equals(c2)){
						score++;
					}
				}
			}
			return ((double)score)/((double)cnt);
		}
	}
	
	public static double norm2vector(double[] v1, double[] v2){
		double dis=0;
		for(int i=0;i<v1.length;i++){
			dis=dis+(v1[i]-v2[i])*(v1[i]-v2[i]);
		}
		return dis;
	}
	
	//v1+v2
	public static double[] vector_add(double[] v1, double[] v2){
		
		double[] v=new double[v1.length];
		for(int i=0;i<v1.length;i++){
			v[i]=v1[i]+v2[i];
		}
		return v;
	}
	
	//v1-v2
	public static double[] vector_minus(double[] v1, double[] v2){
		
		double[] v=new double[v1.length];
		for(int i=0;i<v1.length;i++){
			v[i]=v1[i]-v2[i];
		}
		return v;
	}
	
	//f*v1
	public static double[] vector_factor_mutilple(double[] v1, double f){
		
		double[] v=new double[v1.length];
		for(int i=0;i<v1.length;i++){
			v[i]=v1[i]*f;
		}
		return v;
	}
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		
	}

}
