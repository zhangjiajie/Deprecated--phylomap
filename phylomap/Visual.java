package phylomap;
import java.util.*;
public class Visual {

	/**
	 * @param args
	 */
	HashMap hm=null;
	
	
	public void init(HashMap hm){
		this.hm=hm;
	}
	
	public void setHm(HashMap hm){
		this.hm=hm;
	}
	
	int cnt=1;
	
	public double[] segSpace(double x0, double y0, double x2, double y2, double scale, double tree_dis, double error){
		//this.hm=(HashMap)Name_Coords.clone();
		if(x0>x2){
			double temp0=x0;
			double temp1=y0;
			x0=x2;
			y0=y2;
			x2=temp0;
			y2=temp1;
		}
		
		pnp Pnp=new pnp();
		
		double[] control =new double[3];
		double[] c1=new double[2];
		c1[0]=x0;
		c1[1]=y0;
		double[] c2=new double[2];
		c2[0]=x2;
		c2[1]=y2;
		double Lcord=Pnp.dis(c1, c2);
		double Ltarget=scale*tree_dis;
		double Lstring=Lcord;
		
		double k=-(x2-x0)/(y2-y0);
		double xm=(x0+x2)/2.0;
		double ym=(y0+y2)/2.0;
		
		double step=0.0001;
		double x1=xm;
		double y1=k*(x1-xm)+ym;
		double k2=-1/k;
		
		
		if(k2>=0){
			// assume k2>0
			line line1=new line(k2, y0-k2*x0);
			line line3=new line(k,ym-k*xm);
			double angle1=Math.atan(k2);
			double k3=Math.tan(angle1-Math.PI/4);
			line line2=new line(k3, ym-k3*xm);
			double k4=Math.tan(angle1+Math.PI/4);
			line line4=new line(k4, ym-k4*xm);
			line linea=new line(k4, y2-k4*x2);
			line lineb=new line(k4, y0-k4*x0);
			
			
			if(Lcord>=Ltarget){
				control[0]=x0;
				control[1]=y0;
				control[2]=Lcord;
			}else{
			
				double c_error=Ltarget-Lstring;
				while(c_error>0){
					x1=x1+step;
					y1=k*(x1-xm)+ym;
					Lstring=Pnp.LBezier(x0, y0, x1, y1, x2, y2);
					c_error=Ltarget-Lstring;
				}
				double abs_c_error=Math.sqrt(c_error*c_error);
				while(abs_c_error>error){
					step=step/2;
					if(c_error>0){
						x1=x1+step;
						y1=k*(x1-xm)+ym;
						Lstring=Pnp.LBezier(x0, y0, x1, y1, x2, y2);
						c_error=Ltarget-Lstring;
						
					}else{
						x1=x1-step;
						y1=k*(x1-xm)+ym;
						Lstring=Pnp.LBezier(x0, y0, x1, y1, x2, y2);
						c_error=Ltarget-Lstring;
					}
					abs_c_error=Math.sqrt(c_error*c_error);
				}
				
				double x3=xm-(x1-xm);
				double y3=line3.y(x3);
				line linec=new line(k3, y1-k3*x1);
				line lined=new line(k3, y3-k3*x3);
				
				
				
				int cnt_aera1=0;;
				if(k2>=1){
					cnt_aera1=this.aeraPointNo2(line1, linea, line2, hm);
				}else{
					cnt_aera1=this.aeraPointNo1(line1, linea, line2, hm);
				}
				
				int cnt_aera2=0;
				if(k2>=1){
					cnt_aera2=this.aeraPointNo2(line2, linea, line3, hm);
				}else{
					cnt_aera2=this.aeraPointNo1(line2, linea, line3, hm);
				}
				
				int cnt_aera3=0;
				if(k2>=1){
					cnt_aera3=this.aeraPointNo1(line3, line4, linec, hm);
				}else{
					cnt_aera3=this.aeraPointNo2(line3, line4, linec, hm);
				}
				
				int cnt_aera4=0;
				if(k2>=1){
					cnt_aera4=this.aeraPointNo2(line1, line4, linec, hm);
				}else{
					cnt_aera4=this.aeraPointNo1(line1, line4, linec, hm);
				}
				
				
				int cnt_aera5=0;
				if(k2>=1){
					cnt_aera5=this.aeraPointNo1(line2, line1, lineb, hm);
				}else{
					cnt_aera5=this.aeraPointNo2(line2, lineb, line1, hm);
				}
				
				
				int cnt_aera6=0;
				if(k2>=1){
					cnt_aera6=this.aeraPointNo1(line3, line2, lineb, hm);
				}else{
					cnt_aera6=this.aeraPointNo2(line3, lineb, line2, hm);
				}
				
				int cnt_aera7=0;
				if(k2>=1){
					cnt_aera7=this.aeraPointNo2(lined, line4, line3, hm);
				}else{
					cnt_aera7=this.aeraPointNo1(lined, line4, line3, hm);
				}
				
				int cnt_aera8=0;
				if(k2>=1){
					cnt_aera8=this.aeraPointNo1(lined, line4, line1, hm);
				}else{
					cnt_aera8=this.aeraPointNo2(lined, line4, line1, hm);
				}
				
				
				String direction="A";
				line ldir=line3;
				double ddrect=1.0;
				int aera_sum=cnt_aera2+cnt_aera3;
				
				if((cnt_aera3+cnt_aera4)<aera_sum){
					aera_sum=cnt_aera3+cnt_aera4;
					direction="B";
					ldir=line4;
					if(k2>=1){
						ddrect=1.0;
					}else{
						ddrect=-1.0;
					}
					
				}
				
				if((cnt_aera1+cnt_aera2)<aera_sum){
					aera_sum=cnt_aera1+cnt_aera2;
					direction="C";
					ldir=line2;
					ddrect=1.0;
				}
				
				if((cnt_aera7+cnt_aera8)<aera_sum){
					aera_sum=cnt_aera7+cnt_aera8;
					direction="F";
					ldir=line4;
					if(k2>=1){
						ddrect=-1.0;
					}else{
						ddrect=1.0;
					}
				}
				
				if((cnt_aera6+cnt_aera7)<aera_sum){
					aera_sum=cnt_aera6+cnt_aera7;
					direction="D";
					ldir=line3;
					ddrect=-1.0;
				}
				
				if((cnt_aera6+cnt_aera5)<aera_sum){
					aera_sum=cnt_aera6+cnt_aera5;
					direction="E";
					ldir=line2;
					ddrect=-1.0;
				}
				
				control=ldir.controlPoint(scale*tree_dis, Lcord, xm, error, x0, y0, x2, y2, ddrect);
				
				
				//control[0]=x1;
				//control[1]=y1;
				//control[2]=Lstring;
			}
		}else{
			line line1=new line(k2, y0-k2*x0);
			line line3=new line(k,ym-k*xm);
			double angle1=Math.atan(k2);
			double k3=Math.tan(angle1-Math.PI/4);
			line line2=new line(k3, ym-k3*xm);
			double k4=Math.tan(angle1+Math.PI/4);
			line line4=new line(k4, ym-k4*xm);
			line linea=new line(k4, y2-k4*x2);
			line lineb=new line(k4, y0-k4*x0);
			
			if(Lcord>=Ltarget){
				control[0]=x0;
				control[1]=y0;
				control[2]=Lcord;
			}else{
			
				double c_error=Ltarget-Lstring;
				while(c_error>0){
					x1=x1+step;
					y1=k*(x1-xm)+ym;
					Lstring=Pnp.LBezier(x0, y0, x1, y1, x2, y2);
					c_error=Ltarget-Lstring;
				}
				double abs_c_error=Math.sqrt(c_error*c_error);
				while(abs_c_error>error){
					step=step/2;
					if(c_error>0){
						x1=x1+step;
						y1=k*(x1-xm)+ym;
						Lstring=Pnp.LBezier(x0, y0, x1, y1, x2, y2);
						c_error=Ltarget-Lstring;
						
					}else{
						x1=x1-step;
						y1=k*(x1-xm)+ym;
						Lstring=Pnp.LBezier(x0, y0, x1, y1, x2, y2);
						c_error=Ltarget-Lstring;
					}
					abs_c_error=Math.sqrt(c_error*c_error);
				}
				
				double x3=xm-(x1-xm);
				double y3=line3.y(x3);
				line linec=new line(k3, y1-k3*x1);
				line lined=new line(k3, y3-k3*x3);
				
				
				
				int cnt_aera1=0;
				if(k2<=-1){
					 cnt_aera1=this.aeraPointNo2(line3, linec, line4, hm);
				}else{
					cnt_aera1=this.aeraPointNo1(line3, linec, line4, hm);
				}
				
				int cnt_aera2=0;
				if(k2<=-1){
					cnt_aera2=this.aeraPointNo2(line4, linec, line1, hm);
				}else{
					cnt_aera2=this.aeraPointNo1(line4, linec, line1, hm);
				}
				
				int cnt_aera3=0;
				if(k2<=-1){
					cnt_aera3=this.aeraPointNo1(line1, line2, linea, hm);
				}else{
					cnt_aera3=this.aeraPointNo2(line1, line2, linea, hm);
				}
				
				int cnt_aera4=0;
				if(k2<=-1){
					cnt_aera4=this.aeraPointNo2(line2, line3, linea, hm);
				}else{
					cnt_aera4=this.aeraPointNo1(line3, line2, linea, hm);
				}
				
				int cnt_aera5=0;
				if(k2<=-1){
					cnt_aera5=this.aeraPointNo1(line4, lined, line3, hm);
				}else{
					cnt_aera5=this.aeraPointNo2(line4, lined, line3, hm);
				}
				
				int cnt_aera6=this.aeraPointNo2(line1, lined, line4, hm);
				if(k2<=-1){
					cnt_aera6=this.aeraPointNo1(line1, lined, line4, hm);
				}else{
					cnt_aera6=this.aeraPointNo2(line1, lined, line4, hm);
				}
				
				int cnt_aera7=this.aeraPointNo1(lineb, line1, line2, hm);
				if(k2<=-1){
					cnt_aera7=this.aeraPointNo2(lineb, line2, line1, hm);
				}else{
					cnt_aera7=this.aeraPointNo1(lineb, line2, line1, hm);
				}
				
				int cnt_aera8=this.aeraPointNo1(lineb, line2, line3, hm);
				if(k2<=-1){
					cnt_aera8=this.aeraPointNo1(lineb, line2, line3, hm);
				}else{
					cnt_aera8=this.aeraPointNo2(lineb, line2, line3, hm);
				}
				
				String direction="A";
				line ldir=line3;
				double ddrect=1.0;
				int aera_sum=cnt_aera1+cnt_aera8;
				
				if((cnt_aera1+cnt_aera2)<aera_sum){
					aera_sum=cnt_aera1+cnt_aera2;
					direction="B";
					ldir=line4;
					ddrect=1.0;
				}
				
				if((cnt_aera7+cnt_aera8)<aera_sum){
					aera_sum=cnt_aera7+cnt_aera8;
					direction="C";
					ldir=line2;
					if(k2<=-1){
						ddrect=-1.0;
					}else{
						ddrect=1.0;
					}
				}
				
				if((cnt_aera4+cnt_aera5)<aera_sum){
					aera_sum=cnt_aera4+cnt_aera5;
					direction="F";
					ldir=line3;
					ddrect=-1.0;
				}
				
				if((cnt_aera5+cnt_aera6)<aera_sum){
					aera_sum=cnt_aera5+cnt_aera6;
					direction="D";
					ldir=line4;
					ddrect=-1.0;
				}
				
				if((cnt_aera3+cnt_aera4)<aera_sum){
					aera_sum=cnt_aera3+cnt_aera4;
					direction="E";
					ldir=line2;
					if(k2<=-1){
						ddrect=1.0;
					}else{
						ddrect=-1.0;
					}
				}
				
				control=ldir.controlPoint(scale*tree_dis, Lcord, xm, error, x0, y0, x2, y2, ddrect);
			}
		}
		//System.out.println("oL:"+Lcord+"	tL:"+Ltarget+"	bL:"+Lstring);
		double[] cpp=new double[2];
		cpp[0]=control[0];
		cpp[1]=control[1];
		this.hm.put("control"+cnt++, cpp);
		return control;
	}
	
	int aeraPointNo2(line l1, line l2, line l3, HashMap hm){
		int cnt=0;
		Set s=hm.keySet();
		Iterator it=s.iterator();
		while(it.hasNext()){
			double[] coord=(double[])hm.get(it.next());
			double test1=l1.posTest(coord[0], coord[1]);
			double test2=l2.posTest(coord[0], coord[1]);
			double test3=l3.posTest(coord[0], coord[1]);
			if(test1<=0 && test2<=0 && test3>=0){
				cnt++;
			}
		}
		return cnt;
	}
	
	int aeraPointNo1(line l1, line l2, line l3, HashMap hm){
		int cnt=0;
		Set s=hm.keySet();
		Iterator it=s.iterator();
		while(it.hasNext()){
			double[] coord=(double[])hm.get(it.next());
			double test1=l1.posTest(coord[0], coord[1]);
			double test2=l2.posTest(coord[0], coord[1]);
			double test3=l3.posTest(coord[0], coord[1]);
			if(test1<=0 && test2>=0 && test3>=0){
				cnt++;
			}
		}
		return cnt;
	}
	
	public double[] nonOverlapcoords(double[][] coords, double[] initcoord, int idend){
		double[] nextCoord=new double[2];
		double step=0.0002;
		double a=initcoord[0];
		double b=initcoord[1];
		double m=0;
		double n=0;
		boolean flag=false;
		double s_length=0.2;
		double s_height=0.1;
		
		for(int i=0;i<idend;i++){
			double x=coords[i][0];
			double y=coords[i][1];
			//if find one with in range, break and continue to search
			if(x>=a && x<a+step && y<=b && y>b-step){
				flag=true;
				break;
			}else{
				x=x-s_length;
				y=y-s_height;
				if(x>=a && x<a+step && y<=b && y>b-step){
					flag=true;
					break;
				}	
			}
		}
		
		while(flag){
		
			m=m+1.0;
			n=n+1.0;
			
			
			//-m -n Lower Left
			flag=false;
			for(int i=0;i<idend;i++){
				double x=coords[i][0];
				double y=coords[i][1];
				//if find one with in range, break and continue to search
				if(x<a-(m-1)*step && x>=a-m*step && y<=b-(n-1)*step && y>b-n*step){
					flag=true;
					break;
				}else{
					x=x-s_length;
					y=y-s_height;
					if(x<a-(m-1)*step && x>=a-m*step && y<=b-(n-1)*step && y>b-n*step){
						flag=true;
						break;
					}	
				}
			}
			if(!flag){
				n=-(n-1);
				m=-m;
				
				break;
			}
			
			//-m +n Upper Left 
			flag=false;
			for(int i=0;i<idend;i++){
				double x=coords[i][0];
				double y=coords[i][1];
				//if find one with in range, break and continue to search
				if(x<=a-(m-1)*step && x>=a-m*step && y>=b+(n-1)*step && y<=b+n*step){
					flag=true;
					break;
				}else{
					x=x-s_length;
					y=y-s_height;
					if(x<=a-(m-1)*step && x>=a-m*step && y>=b+(n-1)*step && y<=b+n*step){
						flag=true;
						break;
					}	
				}
			}
			if(!flag){
				m=-m;
			
				break;
			}
			
			//+m +n Upper Right
			flag=false;
			for(int i=0;i<idend;i++){
				double x=coords[i][0];
				double y=coords[i][1];
				//if find one with in range, break and continue to search
				if(x>=a+(m-1)*step && x<=a+m*step && y>=b+(n-1)*step && y<=b+n*step){
					flag=true;
					break;
				}else{
					x=x-s_length;
					y=y-s_height;
					if(x>=a+(m-1)*step && x<=a+m*step && y>=b+(n-1)*step && y<=b+n*step){
						flag=true;
						break;
					}	
				}
			}
			if(!flag){
				
				m=m-1;
				break;
			}
			
			
			//+m -n Lower Right
			flag=false;
			for(int i=0;i<idend;i++){
				double x=coords[i][0];
				double y=coords[i][1];
				//if find one with in range, break and continue to search
				if(x>=a+(m-1)*step && x<=a+m*step && y<=b-(n-1)*step && y>=b-n*step){
					flag=true;
					break;
				}else{
					x=x-s_length;
					y=y-s_height;
					if(x>=a+(m-1)*step && x<=a+m*step && y<=b-(n-1)*step && y>=b-n*step){
						flag=true;
						break;
					}	
				}
			}
			if(!flag){
				n=-n;
				
				break;
			}
			
		}
		
		
		nextCoord[0]=a+m*step;
		nextCoord[1]=b+n*step;
		
		return nextCoord;
	}
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
