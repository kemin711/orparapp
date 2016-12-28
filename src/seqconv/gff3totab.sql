-- ustma1 database created on shake to store all info

drop table if exists source;
create table source (
	genomicid varchar(36) primary key, 
	strand char(1),
	start integer, end integer,
	mol_type varchar(24),
	chromosome integer,
	organism varchar(100) references organism
);

load data local infile '/home/kzhou/work/analysis/externalGenomes/Cryptococcus_neoformans/source.tab'
into table source ignore 1 lines;

create table gene (
	genomicid varchar(36) not null references source, 
	strand char(1), 
	start integer, 
	end integer, 
	id varchar(36) primary key,
	attribute text,
	unique key loc (genomicid,strand,start,end)
);
load data local infile '/home/kzhou/work/analysis/externalGenomes/Cryptococcus_neoformans/gene.tab'
into table gene ignore 1 lines;

drop table if exists transcript;
create table transcript (
	transcriptid varchar(36) primary key,
	RNA_type varchar(12),
	geneid varchar(36) references gene,
	attribute text,
	key gid (geneid)
);
load data local infile '/home/kzhou/work/analysis/externalGenomes/Cryptococcus_neoformans/transcript.tab'
into table transcript ignore 1 lines;

/** different transcripts can share the same exon **/
drop table if exists exon;
create table exon (
	genomicid varchar(36) not null references source, 
	strand char(1), 
	start integer, 
	end integer, 
	exon_number integer, 
	transcriptid varchar(36) references transcript, 
	primary key tid (transcriptid, exon_number),
	key loc (genomicid, strand, start, end)
);

load data local infile '/home/kzhou/work/analysis/externalGenomes/Cryptococcus_neoformans/exon.tab'
into table exon ignore 1 lines;

select * 
from transcript t join gene g on t.geneid=g.id
where t.attribute != g.attribute;
-- only 109 attributes for genes are different from that of transcripts
-- it is redundant because there is not much information from
-- automated annotation programs.


select transcriptid, locus_tag
from exon;
-- 11559

select distinct transcriptid, locus_tag
from exon;
-- 6632

select distinct locus_tag
from exon;
-- 6632

-- it is one-on-one for transcriptid <-> gene

select g.genomicid, g.strand, g.id, g.attribute,
	e.RNA_type, e.transcriptid, e.attributes
from gene g join exon e
	on g.id=e.locus_tag
where g.attribute != e.attributes
;

select transcriptid, count(*)
from exon
group by transcriptid;

-- need to add aa sequence later
-- this version of mysql does not know references syntax!
create table protein (
	proteinid varchar(36) primary key, 
	geneid varchar(36) references gene, 
	product varchar(250), 
	attribute text,
	key gid (geneid)
);
load data local infile '/home/kzhou/work/analysis/externalGenomes/Cryptococcus_neoformans/protein.tab'
into table protein ignore 1 lines;

-- CDS --> protein

drop table if exists CDS;
create table CDS (
	genomicid varchar(36) not null, 
	strand char(1), 
	phase integer, -- 0, 1, and 2
	start integer, 
	end integer, 
	id varchar(70), -- all mRNA has NULL id!
	exon_number integer, 
	proteinid varchar(36) references protein, 
	primary key pid (proteinid, exon_number),
	key loc (genomicid, strand, start, end)
);
load data local infile '/home/kzhou/work/analysis/externalGenomes/Cryptococcus_neoformans/CDS.tab'
into table CDS ignore 1 lines;

-- start codon

/** codon is always 3 bases,
	But there are cases where the codong was split
	into multiple exons, in this case there might
	be two start codong locations:

	Example: protein XP_566650.1
CDS: on - strand:
	8954-9010 9057-9181 9240-9714 9770-9851 9904-9905

	The start codon is 
	start_codon 679904  679905  .   -   0   exnum=1 
	start_codon 679851  679851  .   -   1   exon_number=2

**/
drop table if exists start_codon;
create table start_codon (
	genomicid varchar(36), 
	strand char(1), 
	start integer, end integer, 
	proteinid varchar(36) references protein,
	primary key (proteinid, start),
	key loc (genomicid,strand,start),
	key pid (proteinid)
);
load data local infile '/home/kzhou/work/analysis/externalGenomes/Cryptococcus_neoformans/start_codon.tab'
into table start_codon ignore 1 lines;

-- not unique
select genomicid, strand,start, count(*) as count
from start_codon
group by genomicid, strand,start
having count(*)>1;

/** how could the same protein have two different start
codons?
**/

select *
from (
	select proteinid, count(*) as count
	from start_codon
	group by proteinid having count(*) > 1
	) a join start_codon s on a.proteinid=s.proteinid
;

create table stop_codon (
	genomicid varchar(36), 
	strand char(1), 
	start integer, end integer, 
	proteinid varchar(36) references protein,
	primary key (proteinid, start),
	key pid (proteinid),
	key loc (genomicid,strand,start)
);

load data local infile '/home/kzhou/work/analysis/externalGenomes/Cryptococcus_neoformans/stop_codon.tab'
into table stop_codon ignore 1 lines;

select genomicid, strand, start, count(*)
from stop_codon
group by genomicid, strand, start
having count(*)>1;



-- which proteins don't have start codons?
select count(*)
from protein p join start_codon s on p.proteinid=s.proteinid
;

select p.*
from protein p left outer join start_codon s on p.proteinid=s.proteinid
where s.proteinid is null
;

-- those without stop codons
select p.*
from protein p left outer join stop_codon s on p.proteinid=s.proteinid
where s.proteinid is null
;
-- all got stop codons

-- the source maps the chromosome location of each DNA fragment
-- not sure where is the assembly instruction, for gene based 
-- analysis this is not important.

create table organism (
	name varchar(100) primary key,
	taxon integer unsigned,
	attribute text
);
insert into organism values 
('Cryptococcus neoformans var. neoformans JEC21',
 214684,
 'strain=JEC21; variety=neoformans; serotype=D; note=MAT-alpha');

/*
insert into organism values 
('Ustilago maydis', 'db_xref=taxon:5270; organelle=mitochondrion'),
('Ustilago maydis 521', 'db_xref=taxon:237631; strain=521');
*/


-- now wee need to load the protein sequence into the database
-- we can use the fasta for blast directly from the original 
-- download

-- protein sequences

drop table proteinseq;
create table proteinseq (
	gi integer primary key,
	proteinid varchar(36) unique not null,
	title text,
	sequence longtext not null
);

load data local infile 
'/home/kzhou/work/analysis/externalGenomes/Cryptococcus_neoformans/proteinseq.tab'
--into table proteinseq ignore 1 lines;
into table proteinseq;

-- the genomicseq table could be combined with source.
-- it is fine to have some redundancy and isolation.
drop table if exists genomicseq;
create table genomicseq (
	gi integer primary key,
	genomicid varchar(36) unique not null,
	title text,
	sequence longtext not null
);

load data local infile 
'/home/kzhou/work/analysis/externalGenomes/Cryptococcus_neoformans/genomicseq.tab'
-- into table genomicseq ignore 1 lines;
into table genomicseq;

--- protein table transformation

alter  table protein rename to protein_raw;

drop table if exists protein;
create table protein (
	id integer primary key,
	accesssion varchar(36) unique,
	geneid varchar(36) references gene,
	product varchar(250),
	attribute text,
	seq longtext
);

insert into protein
select s.gi, p.proteinid, p.geneid, p.product, p.attribute,
	s.seq
from protein_raw p join proteinseq s on p.proteinid=s.proteinid
;
