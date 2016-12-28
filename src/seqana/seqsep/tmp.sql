show index from animalselfblp_need_dynamic;
create table animal_nd (
	query integer unsigned,
	target integer unsigned,
	primary key(query,target),
	key query(query),
	key target(target)
);
insert into animal_nd
select *
from animalselfblp_need_dynamic;

