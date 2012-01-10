package phylomap;

public class line {

	/**
	 * @param args
	 */
	double k=1;
	double b=0;
	public line(double k, double  b){
		this.k=k;
		this.b=b;
	}
	
	public double posTest(double x, double y){
		double test=0;
		
		test=y-k*x-b; 
		
		return test;
	}
	
	public double y(double x){
		double y=0;
		y=k*x+b;
		
		return y;
		
	}
	
	public double[] controlPoint(double Ltarget, double Lcoord, double x1, double error, double x0, double y0, double x2, double y2, double direction){
		double step=0.0001*direction;
		//double step=0.1;
		double y1=0;
		pnp Pnp=new pnp();
		double Lstring=Lcoord;
		double[] cp=new double[3];
		double c_error=Ltarget-Lstring;
		while(c_error>0){
			x1=x1+step;
			y1=this.y(x1);
			Lstring=Pnp.LBezier(x0, y0, x1, y1, x2, y2);
			c_error=Ltarget-Lstring;
		}
		double abs_c_error=Math.sqrt(c_error*c_error);
		while(abs_c_error>error){
			step=step/2;
			if(c_error>0){
				x1=x1+step;
				y1=this.y(x1);
				Lstring=Pnp.LBezier(x0, y0, x1, y1, x2, y2);
				c_error=Ltarget-Lstring;
				
			}else{
				x1=x1-step;
				y1=this.y(x1);
				Lstring=Pnp.LBezier(x0, y0, x1, y1, x2, y2);
				c_error=Ltarget-Lstring;
			}
			abs_c_error=Math.sqrt(c_error*c_error);
		}
		
		cp[0]=x1;
		cp[1]=y1;
		cp[2]=Lstring;
		//System.out.println(abs_c_error);
		return cp;
		
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		System.out.println("Hellow world in linux!");
	}

}
