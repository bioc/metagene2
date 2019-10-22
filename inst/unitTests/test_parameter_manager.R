
test_ph = metagene2:::parameter_manager$new(
    param_values=list(A="a", B="a", C="c"),
    param_validations=list(
        A=function(x) stopifnot(x=="a"),
        B=function(x) stopifnot(x %in% c("a", "b"))),
    overall_validation=function(params) stopifnot(params$A==params$B))

test.parameter_manager_get_set = function() {
    test = test_ph$clone()
    
    # Test get, set.
    checkTrue(test$get("C")=="c")
    test$set("C", "foo")
    checkTrue(test$get("C")=="foo")
}

test.parameter_manager_update_params = function() {
    test = test_ph$clone()

    # Test update_params
    test$update_params(C="bar")
    checkTrue(test$get("C")=="bar")
    
    # Test update params with name inference
    C=42
    test$update_params(C)
    checkTrue(test$get("C")==42)
}

test.parameter_manager_have_params_changed = function() {
    test = test_ph$clone()

    # Test have_params_changed TRUE/FALSE
    checkTrue(test$have_params_changed(C="d"))
    checkTrue(!test$have_params_changed(C="c"))
    
    # Test have_params_changed with name inference
    C="d"
    checkTrue(test$have_params_changed(C))
    C="c"
    checkTrue(!test$have_params_changed(C))
    
}

test.parameter_manager_rollback = function() {
    test = test_ph$clone()
        
    # Make sure rollback works when single parameter validation fails.
    old_B = test$get("B")
    
    obs <- tryCatch(test$set("B", "Patate"), error = conditionMessage)
    checkTrue(obs == 'x %in% c("a", "b") is not TRUE')
    checkTrue(test$get("B") == old_B)
    
    # Make sure rollback works when combined validation fails (with set)
    obs <- tryCatch(test$set("B", "b"), error = conditionMessage)
    checkTrue(obs == 'params$A == params$B is not TRUE')
    checkTrue(test$get("B") == old_B)    
    
    # Make sure rollback works when combined validation fails (with update_params)
    old_C = test$get("C")    
    obs <- tryCatch(test$update_params(B="b", C="WRONG"), error = conditionMessage)
    checkTrue(obs == 'params$A == params$B is not TRUE')
    checkTrue(test$get("B") == old_B)        
    checkTrue(test$get("C") == old_C)        
}

test.parameter_manager_lock = function() {
    test = metagene2:::parameter_manager$new(
        param_values=list(A="a"),
        locked=TRUE)
        
    obs <- tryCatch(test$set("G", 7), error = conditionMessage)
    checkTrue(obs == "Trying to access unknown parameter G")

    test = metagene2:::parameter_manager$new(
        param_values=list(A="a"),
        locked=FALSE)
        
    test$set("G", 7)
    checkTrue(test$get("G")==7)
}
