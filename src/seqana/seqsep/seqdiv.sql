/* trying to divide the nr into animals
*/
select id from tax_division where name in('Invertebrates', 'Mammals', 'Primates', 'Rodents', 'Vertebrates');
+----+
| id |
+----+
|  1 |
|  2 |
|  5 |
|  6 |
| 10 |
+----+


select *
from tax_nodes
where division in (1,2,5,6,10);

-- those gi for all animals
select distinct x.gi
from prtid2taxid x join tax_nodes n on x.taxid=n.id
where n.division in (1,2,5,6,10)
;

-- all proteins from animal
drop view animalgi;

create table animalgi as
select distinct x.gi
from prtgi2taxid x join tax_nodes n on x.taxid=n.id
where n.division in (1,2,5,6,10)
;
alter table animalgi add primary key(gi);

optimize table nr;

create table animalnr (
	id integer unsigned primary key,
	title text,
	length integer,
	sequence longtext
);
insert into animalnr
select nr.id, nr.title, nr.length, nr.sequence
from animalgi g join nr on g.gi=nr.id
;

/** now plants */
select * from tax_division;


create table plantgi as
select distinct x.gi
from prtgi2taxid x join tax_nodes n on x.taxid=n.id
where n.division=4;

alter table plantgi add primary key(gi);

create table plantnr (
	id integer unsigned primary key,
	title text,
	length integer,
	sequence longtext
);

insert into plantnr
select nr.id, nr.title, nr.length, nr.sequence
from plantgi p join nr on p.gi=nr.id;


