/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package processpileup;

/**
 *
 * @author Radoslaw Suchecki radoslaw.suchecki@adelaide.edu.au
 */
public class BaseCall implements Comparable<BaseCall> {

    private final String base;
    private final int count;

    public BaseCall(String base, int count) {
        this.base = base;
        this.count = count;
    }

    public String getBase() {
        return base;
    }

    public int getCount() {
        return count;
    }

    @Override
    public int compareTo(BaseCall another) {
        return another.getCount() - getCount();
    }

}
