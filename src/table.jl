# If we add Tables.jl as a dependency, we can simplify some code here and in params.jl.

function _to_named_tuple(table)
	names = propertynames(table)
	values = (getproperty(table,name) for name in names)
	NamedTuple{(names...,)}(values)
end


function subset_rows(table, rows)
	names = propertynames(table)
	values = (getproperty(table,name)[rows] for name in names)
	NamedTuple{(names...,)}(values)
end

remove_columns(table::NamedTuple, cols) = Base.structdiff(table,NamedTuple{cols})
remove_columns(table, cols) = remove_columns(_to_named_tuple(table), cols)

hcat_tables(t1::NamedTuple, t2::NamedTuple) = merge(t1,t2)
hcat_tables(t1::NamedTuple, t2) = merge(t1,_to_named_tuple(t2))
hcat_tables(t1, t2) = merge(_to_named_tuple(t1),t2)
