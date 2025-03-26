
%................................................................... PRIVATE
%.........................................................................

private define get_params ()
{
   variable i, n = get_num_pars ();
   variable p = Struct_Type [n];
   
   for (i = 1; i <= n; i++)
     {
    p[i-1] = get_par_info (i);

     }
   
   return p;
}

public define how_many_gaussian ()
{
    variable n, info_par, last_gaussian, num;
    n = get_num_pars();
    info_par = get_par_info(n);
    last_gaussian = info_par.name;

    sscanf(last_gaussian, "egauss(%d).area", &num);

    return num;
}


private define test_params (particle, evalfun)
{
    variable info, nbins, free_pars, chisquared;
    set_params(particle);         % Apply particle parameters to the system
    () = @evalfun(&info);         % Evaluate; evalfun fills in 'info'
    
    nbins = info.num_bins;
    free_pars = info.num_variable_params;
    chisquared = info.statistic / (nbins - free_pars);
    return chisquared;            % Return reduced chi-square
}

%...................................................................PUBLIC
%.........................................................................
public define create_line_model(name) {
    write_plot("spec");
    variable cmd;
    cmd = "python " + MODEL_DIR + "/lines.py " + name;
    system(cmd);

    cmd = "python " + MODEL_DIR + "/lines.py";
    system(cmd);
    
    () = evalfile("set_line_model.sl");
    vmessage("Lines saved in model under name linemodel");
}

    
public define fit_line_model(evalfun) {

    () = evalfile("set_line_parameters.sl");

    variable num = how_many_gaussian();
    
    if (num <= 0) {
       vmessage("Invalid number of gaussians: %i", num);
       return;
    }
    
    variable best_p = get_params();
    variable best_stat = 1.0e99;

    variable i;
    vmessage("number of line candidates detected: %i", num);
    
    for (i = 1; i < num+1; i++) {
        vmessage("iteration number %f of %f",i,num);
    
        variable p = get_params();
        variable stat = test_params(p, evalfun);

        variable paramName_area = sprintf("egauss(%i).area", i);
        variable paramName_sigma = sprintf("egauss(%i).sigma", i);
        variable paramName_center = sprintf("egauss(%i).center", i);
        
        thaw(paramName_area);
        thaw(paramName_sigma);
        thaw(paramName_center);
        
        fit_counts();
        variable new_p = get_params();
        variable new_stat = test_params(new_p, evalfun);
        
        
        vmessage("next stat: %f", new_stat);
        
        if (new_stat < stat-0.005) {
            vmessage("line accepted");

            best_stat = new_stat;
            best_p = new_p;
            set_params(new_p);
            }
        else{
            set_par(paramName_area, 0);
            freeze(paramName_area);
            freeze(paramName_sigma);
            freeze(paramName_center);
            }
    }
    

  set_params(best_p);

}
