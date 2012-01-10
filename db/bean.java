package db;


public class bean {

	/**
	 * @param args
	 */
	String id="";
	String sequence="";
	String annotation="";
	public double x;
	public double y;
	public double[] acc_vector;
	
	public void setId(String id){
		this.id=id;
	}
	public String getId(){
		return this.id;
	}
	
	public void setSequence(String sequence){
		this.sequence=sequence;
	}
	public String getSequence(){
		return this.sequence;
	}
	
	public void setAnnotation(String annotation){
		this.annotation=annotation;
	}
	public String getAnnotation(){
		return this.annotation;
	}
	
	public String getCoords(){
		return this.x+"		"+this.y;
	}
	
	public String toString(){
		return this.id+"	"+this.annotation+"	"+this.sequence+"\n";
	}
	public boolean equals(Object o){
		/*
		if(!o.getClass().isInstance(this)){
			return false;
		}*/
		if(o instanceof bean){
			bean O=(bean)o;
			if(this==O){
				return true;
			}else{
				if(O.getId().equals(this.id)){
					return true;
				}else{
					return false;
				}
			}
		}else{
			return false;
		}
	
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
