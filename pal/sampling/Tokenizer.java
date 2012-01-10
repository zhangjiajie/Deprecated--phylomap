package pal.sampling;
import java.io.*;
import java.util.*;

public class Tokenizer {
	private String input;
	private LinkedList<String> tlist=new LinkedList<String>();
	private int index=0;
	public Tokenizer(String inputs){
		this.input=input;
		String[] ss=inputs.split("\\s");
		for(int i=0; i<ss.length; i++){
			if(!"".equals(ss[i].trim())){
				tlist.add(ss[i]);
			}
		}
	}
	public String nextToken(){
		if(index<tlist.size()){
			return tlist.get(index++);
		}else{
			return null;
		}
	}
	
	public static void main(String[] args) {
		String s="t2         0.000000 0.453449 0.471797 0.442784 0.532411 0.569423 0.569717";
		Tokenizer tkzer=new Tokenizer(s);
		String tk=tkzer.nextToken();
		while(tk!=null){
			
			System.out.println(tk);
			tk=tkzer.nextToken();
		}
	}
	
}
