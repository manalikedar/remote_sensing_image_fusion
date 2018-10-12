package com.cg.project.bean;

import com.cg.project.exception.InvalidNuumberRangeException;

public interface MathServices {
	public abstract int addNums(int n1,int n2) throws InvalidNuumberRangeException;
	abstract int subNums(int n1,int n2) throws InvalidNuumberRangeException;
	int multiNums(int n1,int n2) throws InvalidNuumberRangeException;
}
