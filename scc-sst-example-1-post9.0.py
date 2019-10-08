import sys
import sst
import ConfigParser, argparse

#####################################################################
# This example is designed to work with SST releases AFTER 9.0
# Namely, the github master branch and the upcoming 9.1 release
# It will produce errors if run with the 9.0 release due to changes
# in how subcomponents are specified
#####################################################################

# Command line args
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--statfile", help="stat output file", default="./example1.csv")
parser.add_argument("-c", "--clock", help="core frequency", default="2.5GHz")
parser.add_argument("-l3", "--l3size", help="size of l3, with units, SI units ok. Ki and K both treated as ^2.", default="1MiB")
parser.add_argument("-s", "--singlestream", help="number of cores running the singlestream generator", default="20")
args = parser.parse_args();

# Read command line args
statFile = args.statfile
core_clock = args.clock
l3_size = args.l3size
singlestreamcpus = int(args.singlestream)

#quiet = True
quiet = False

# Debug parameters
debugAll = 0
debugL1 = max(debugAll, 0)
debugL2 = max(debugAll, 0)
debugL3 = max(debugAll, 0)
debugMemCtrl = max(debugAll, 0)
debugNIC = max(debugAll, 0)
debugLev = 3
memHVerbose = 1 # 0 = no output, 1 = print warnings, 2 = warnings and dump state on emergency shutdown (fatal/kill/etc.)

################## Problem size configuration #########################
# Some cores run singlestream, the rest of the cores run SpMV

corecount = 28

# Little bit of error checking
if singlestreamcpus > corecount:
    print "Error: corecount is only " + str(corecount) + " but config asked for " + str(singlestreamcpus) + " cores to run SingleStreamGenerator"
    sys.exit(0)


# SingleStream parameters
ss_total_elems = 5000000    # Total number of elements to stream among all cores
ss_length = 16              # Size of each element in bytes

# SPMV parameters
nx_per_core = 2000  # Cols per core
ny_per_core = 2000  # Rows per core
nnz = 12            # Non-zeros per row


################## Architecture configuration #########################
# Memory
memCapacity = 8     # In GB
memPageSize = 4     # In KB
memNumPages = memCapacity * 1024 * 1024 / memPageSize
coherence_protocol = "MESI"

# NoC
mesh_stops_x    = 6
mesh_stops_y    = 5
mesh_clock      = 1500   #MHz
ctrl_mesh_flit      = 8
data_mesh_flit      = 36
mesh_link_latency   = "150ps"
ctrl_mesh_link_bw   = str( mesh_clock * 1000 * 1000 * ctrl_mesh_flit ) + "B/s"
data_mesh_link_bw   = str( mesh_clock * 1000 * 1000 * data_mesh_flit ) + "B/s"
ctrl_network_buffers = "32B"
data_network_buffers = "288B"

# NoC/MemNIC parameters
ctrl_network_params = {
        "link_bw" : ctrl_mesh_link_bw,
        "flit_size" : str(ctrl_mesh_flit) + "B",
        "input_buf_size" : ctrl_network_buffers,
        "route_y_first" : 1,
}

data_network_params = {
        "link_bw" : data_mesh_link_bw,
        "flit_size" : str(data_mesh_flit) + "B",
        "input_buf_size" : data_network_buffers,
        "route_y_first" : 1
}

nic_params = {
    "accept_region" : 0,
    "debug" : debugNIC,
    "debug_level" : debugLev,
}

ctrl_nic_params = {
    "link_bw" : ctrl_mesh_link_bw,
    "in_buf_size" : ctrl_network_buffers,
    "out_buf_size" : ctrl_network_buffers,
}

data_nic_params = {
    "link_bw" : data_mesh_link_bw,
    "in_buf_size" : data_network_buffers,
    "out_buf_size" : data_network_buffers,
}

# L1 cache parameters
l1_cache_params = {
    "cache_frequency"    : core_clock,
    "coherence_protocol" : coherence_protocol,
    "replacement_policy" : "lru",
    "cache_size"         : "32KiB",
    "associativity"      : 8,
    "cache_line_size"    : 64,
    "access_latency_cycles" : 3,
    "tag_access_latency_cycles" : 1,
    "mshr_num_entries"   : 12,
    "maxRequestDelay"    : 1000000,
    "events_up_per_cycle" : 2,
    "mshr_latency_cycles" : 2,
    "max_requests_per_cycle" : 4,
    "L1"                 : 1,
    "debug"              : debugL1,
    "debug_level"        : debugLev,
    "verbose"            : memHVerbose,
}

# L2 prefetcher parameters
l2_prefetch_params = {
    "reach" : 16,
    "detect_range" : 1
}

# L2 cache parameters
l2_cache_params = {
    "cache_frequency"    : core_clock,
    "coherence_protocol" : coherence_protocol,
    "replacement_policy" : "lru",
    "cache_size"         : "512KiB",
    "associativity"      : 16,
    "cache_line_size"    : 64,
    "access_latency_cycles" : 8,
    "tag_access_latency_cycles" : 3,
    "mshr_num_entries"   : 16,
    "mshr_latency_cycles" : 3,
    "response_link_width" : "72B",
    "debug"              : debugL2,
    "debug_level"        : debugLev,
    "verbose"            : memHVerbose,
}

# L3 cache parameters
l3_cache_params = {
    "cache_frequency" : str(mesh_clock) + "MHz",
    "coherence_protocol" : coherence_protocol,
    "replacement_policy" : "random",
    "cache_size" : l3_size,
    "associativity" : 16,
    "cache_line_size" : 64,
    "access_latency_cycles" : 20,
    "tag_latency_cycles" : 5,
    "mshr_num_entries" : 32,
    "mshr_latency_cycles" : 5,
    "num_cache_slices" : corecount,
    "slice_allocation_policy" : "rr",
    "cache_type" : "noninclusive_with_directory",
    "noninclusive_directory_entries" : 65536,
    "noninclusive_directory_associativity" : 16,
    "debug" : debugL3,
    "debug_level" : debugLev,
    "verbose"       : memHVerbose,
}

# MemoryController parameters
memctrl_params = {
    "clock"     : "1.2GHz",
    "verbose"   : memHVerbose,
    "max_requests_per_cycle" : 3,
    "backing" : "none",
    "cache_line_size" : 64,
    "debug" : debugMemCtrl,
    "debug_level" : 3,
}

# MemoryConroller backend paramters
membackend_params = {
    "id" : 0,
    "addrMapper" : "memHierarchy.roundRobinAddrMapper",
    "channels" : 2,
    "channel.transaction_Q_size" : 32,
    "channel.numRanks" : 2,
    "channel.rank.numBanks" : 16,
    "channel.rank.bank.CL" : 15,
    "channel.rank.bank.CL_WR" : 12,
    "channel.rank.bank.RCD" : 15,
    "channel.rank.bank.TRP" : 15,
    "channel.rank.bank.dataCycles" : 4,
    "channel.rank.bank.pagePolicy" : "memHierarchy.simplePagePolicy",
    "channel.rank.bank.transactionQ" : "memHierarchy.reorderTransactionQ",
    "channel.rank.bank.pagePolicy.close" : 0,
    "printconfig" : 0,
    "channel.printconfig" : 0,
    "channel.rank.printconfig" : 0,
    "channel.rank.bank.printconfig" : 0,
}

# CPU parameters
miranda_params = {
    "verbose" : 0,
    "clock" : core_clock,
    "max_reqs_cycle" : 2,
    "max_reorder_lookups" : 32,
    "maxmemreqpending" : 20,
    "pagecount" : memNumPages,
    "pagesize" : memPageSize * 1024
}

# SingleStream generator parameters
ss_core_elems = int(ss_total_elems / singlestreamcpus)  # Number of elements to touch
ss_end = ss_core_elems * ss_length * singlestreamcpus # Last address accessed by the singlestream cpus
singlestream_gen_params = {
    "verbose" : 0,
    "count" : ss_core_elems,
    "length" : ss_length,
    "max_address" : memCapacity * 1024 * 1024 * 1024 - 16 # Max address in system
}

print "Configuring " + str(singlestreamcpus) + " cores to run singlestream"
print " - Number of elements accessed per core: " + str(ss_core_elems)
print " - Addresses touched: 0 - " + str(ss_end)

# SpMV generator parameters
matrix_nx = nx_per_core * (corecount - singlestreamcpus)
matrix_ny = ny_per_core * (corecount - singlestreamcpus)
lhs_start = ss_end + (ss_end % 64) # Pad so we start on a cache line boundary
rhs_start = lhs_start + 8 * matrix_ny
row_start = rhs_start + 8 * matrix_ny
col_start = row_start + 8 * (matrix_nx + 1)
elements_start = col_start + 8 * (nnz * matrix_nx)
spmv_gen_params = {
    "verbose" : 1,
    "matrix_nx" : matrix_nx,
    "matrix_ny" : matrix_ny,
    "element_width" : 8,
    "lhs_start_addr" : lhs_start,
    "rhs_start_addr" : rhs_start,
    "ordinal_width" : 8,
    "matrix_row_indices_start_addr" : row_start,
    "matrix_col_indices_start_addr" : col_start,
    "matrix_element_start_addr" : elements_start,
    "iterations" : 3,
    "matrix_nnz_per_row" : nnz,
}

print ""
print "Configuring " + str(corecount - singlestreamcpus) + " cores to run SpMV"
print " - nx x ny: " + str(matrix_nx) + " x " + str(matrix_ny)
print " - Rows per core: " + str(ny_per_core)
print " - NonZero/row: " + str(nnz)
print " - LHS vector address: " + str(lhs_start)
print " - RHS vector address: " + str(rhs_start)
print " - Row vector address: " + str(row_start)
print " - Col vector address: " + str(col_start)
print " - Element vector address: " + str(elements_start)


# Build the memories
class MemBuilder:
    def __init__(self, memCapacity, startAddr, memCount):
        self.next_mem_id = 0
        self.mem_count = memCount
        self.mem_capacity = memCapacity
        self.start_addr = startAddr

    def build(self, nodeID):
        if not quiet:
            print "Creating Mem controller " + str(self.next_mem_id) + " out of " + str(self.mem_count) + " on node " + str(nodeID) + "..."
            print " - Capacity: " + str(self.mem_capacity / self.mem_count) + " per Memory."

        mem = sst.Component("mem_" + str(self.next_mem_id), "memHierarchy.MemController")
        mem.addParams(memctrl_params)
        
        membk = mem.setSubComponent("backend", "memHierarchy.timingDRAM")
        membk.addParams(membackend_params)
        membk.addParams({ "mem_size" : str(self.mem_capacity / self.mem_count) + "B" })

        memLink = sst.Link("mem_link_" + str(self.next_mem_id))

        # Define memory NIC
        mlink = mem.setSubComponent("cpulink", "memHierarchy.MemNICFour")
        data = mlink.setSubComponent("data", "kingsley.linkcontrol")
        req = mlink.setSubComponent("req", "kingsley.linkcontrol")
        fwd = mlink.setSubComponent("fwd", "kingsley.linkcontrol")
        ack = mlink.setSubComponent("ack", "kingsley.linkcontrol")
        mlink.addParams(nic_params)
        mlink.addParams({"group" : 3})
        data.addParams(data_nic_params)
        req.addParams(ctrl_nic_params)
        fwd.addParams(ctrl_nic_params)
        ack.addParams(ctrl_nic_params)

        self.next_mem_id = self.next_mem_id + 1
        return (req,"rtr_port", mesh_link_latency), (ack, "rtr_port", mesh_link_latency), (fwd, "rtr_port", mesh_link_latency), (data, "rtr_port", mesh_link_latency)

# Build the L3 slices, one per chip
class L3Builder:
    def __init__(self):
        self.next_l3_id = 0

    def build(self):

        ## L3 ##
        l3Cache = sst.Component("l3cache_" + str(self.next_l3_id), "memHierarchy.Cache")
        l3Cache.addParams(l3_cache_params)
        l3Cache.addParams({ "slice_id" : self.next_l3_id })

        l3nic = l3Cache.setSubComponent("cpulink", "memHierarchy.MemNICFour")
        l3nic.addParams(nic_params)
        l3nic.addParams({"group" : 2})
        l3data = l3nic.setSubComponent("data", "kingsley.linkcontrol")
        l3req = l3nic.setSubComponent("req", "kingsley.linkcontrol")
        l3fwd = l3nic.setSubComponent("fwd", "kingsley.linkcontrol")
        l3ack = l3nic.setSubComponent("ack", "kingsley.linkcontrol")
        l3data.addParams(data_nic_params)
        l3req.addParams(ctrl_nic_params)
        l3fwd.addParams(ctrl_nic_params)
        l3ack.addParams(ctrl_nic_params)

        self.next_l3_id = self.next_l3_id + 1

        return (l3req, "rtr_port", mesh_link_latency), (l3ack, "rtr_port", mesh_link_latency), (l3fwd, "rtr_port", mesh_link_latency), (l3data, "rtr_port", mesh_link_latency)


# Build the cores, 28 cores
class CoreBuilder:
    def __init__(self):
        self.next_core_id = 0

    def build(self):
        ## Core ##
        core = sst.Component("cpu." + str(self.next_core_id), "miranda.BaseCPU")
        core.addParams(miranda_params)

        # Not all cores are running the same thing
        if self.next_core_id < singlestreamcpus:
            gen = core.setSubComponent("generator", "miranda.SingleStreamGenerator")
            gen.addParams(singlestream_gen_params)
            gen.addParams( { "startat" : ss_core_elems * 16 * self.next_core_id } )
            print "  Configuring core " + str(self.next_core_id) + " for SingleStream"
        else:
            gen = core.setSubComponent("generator", "miranda.SPMVGenerator")
            gen.addParams(spmv_gen_params)
            gen.addParams({
                "local_row_start" : ny_per_core * (self.next_core_id - singlestreamcpus),
                "local_row_end" : ny_per_core * (self.next_core_id + 1 - singlestreamcpus) - 1
            })
            print "  Configuring core " + str(self.next_core_id) + " for SpMV"

        ## L1 ##
        l1Cache = sst.Component("l1cache_" + str(self.next_core_id), "memHierarchy.Cache")
        l1Cache.addParams(l1_cache_params)

        ## L2 ##
        l2Cache = sst.Component("l2cache_" + str(self.next_core_id), "memHierarchy.Cache")
        l2Cache.addParams(l2_cache_params)
        ## -> Prefetcher
        pref = l2Cache.setSubComponent("prefetcher", "cassini.StridePrefetcher")
        pref.addParams(l2_prefetch_params)
        ## -> Links
        l2cpu = l2Cache.setSubComponent("cpulink", "memHierarchy.MemLink")
        l2mem = l2Cache.setSubComponent("memlink", "memHierarchy.MemNICFour")
        l2mem.addParams(nic_params)
        l2mem.addParams({"group" : 1})
        l2data = l2mem.setSubComponent("data", "kingsley.linkcontrol")
        l2req = l2mem.setSubComponent("req", "kingsley.linkcontrol")
        l2fwd = l2mem.setSubComponent("fwd", "kingsley.linkcontrol")
        l2ack = l2mem.setSubComponent("ack", "kingsley.linkcontrol")
        l2data.addParams(data_nic_params)
        l2req.addParams(ctrl_nic_params)
        l2fwd.addParams(ctrl_nic_params)
        l2ack.addParams(ctrl_nic_params)

        # Connect L1 & L2
        l1Tol2 = sst.Link("l1tol2_" + str(self.next_core_id))
        l1Tol2.connect( (l1Cache, "low_network_0", "100ps"), (l2cpu, "port", "100ps") )
        l1Tol2.setNoCut()

        # Connect Core & l1
        cpuTol1 = sst.Link("cpulink_" + str(self.next_core_id))
        cpuTol1.connect( (core, "cache_link", "75ps"), (l1Cache, "high_network_0", "75ps") )
        cpuTol1.setNoCut()

        self.next_core_id = self.next_core_id + 1

        return (l2req, "rtr_port", mesh_link_latency), (l2ack, "rtr_port", mesh_link_latency), (l2fwd, "rtr_port", mesh_link_latency), (l2data, "rtr_port", mesh_link_latency)

# Builder instances
coreBuilder = CoreBuilder()
l3Builder = L3Builder()
memBuilder  = MemBuilder(memCapacity * 1024 * 1024 * 1024, 0, 2) # 2 memories

#######################################################
# Build Config
#######################################################
# 28 cores and 30 mesh stops (6x5)
#   Mem = memory
#   Core = CPU + L1 + L2 + L3 slice
#
#   Core -- Core -- Core -- Core -- Core -- Core
#    |       |       |       |       |       | 
#   Core -- Core -- Core -- Core -- Core -- Core
#    |       |       |       |       |       | 
#   Mem  -- Core -- Core -- Core -- Core -- Mem 
#    |       |       |       |       |       | 
#   Core -- Core -- Core -- Core -- Core -- Core
#    |       |       |       |       |       | 
#   Core -- Core -- Core -- Core -- Core -- Core
#
#######################################################

def setNode(nodeId, rtrreq, rtrack, rtrfwd, rtrdata):

    if nodeId == 12 or nodeId == 17:
        memreq, memack, memfwd, memdata = memBuilder.build(nodeId)
        reqport0 = sst.Link("krtr_req_port0_" + str(nodeId))
        fwdport0 = sst.Link("krtr_fwd_port0_" + str(nodeId))
        ackport0 = sst.Link("krtr_ack_port0_" + str(nodeId))
        dataport0 = sst.Link("krtr_data_port0_" + str(nodeId))
        reqport0.connect( (rtrreq, "local0", mesh_link_latency), memreq)
        fwdport0.connect( (rtrfwd, "local0", mesh_link_latency), memfwd)
        ackport0.connect( (rtrack, "local0", mesh_link_latency), memack)
        dataport0.connect( (rtrdata, "local0", mesh_link_latency), memdata)
    else:
        # Place cores on all other routers
        corereq, coreack, corefwd, coredata = coreBuilder.build()
        reqport0 = sst.Link("krtr_req_port0_" + str(nodeId))
        reqport0.connect( (rtrreq, "local0", mesh_link_latency), corereq )
        ackport0 = sst.Link("krtr_ack_port0_" + str(nodeId))
        ackport0.connect( (rtrack, "local0", mesh_link_latency), coreack )
        fwdport0 = sst.Link("krtr_fwd_port0_" + str(nodeId))
        fwdport0.connect( (rtrfwd, "local0", mesh_link_latency), corefwd )
        dataport0 = sst.Link("kRtr_data_port0_" + str(nodeId))
        dataport0.connect( (rtrdata, "local0", mesh_link_latency), coredata )

        # Place L3 slices on all routers
        l3req, l3ack, l3fwd, l3data = l3Builder.build()
        reqport1 = sst.Link("krtr_req_port1_" + str(nodeId))
        reqport1.connect( (rtrreq, "local1", mesh_link_latency), l3req)
        ackport1 = sst.Link("krtr_ack_port1_" + str(nodeId))
        ackport1.connect( (rtrack, "local1", mesh_link_latency), l3ack)
        fwdport1 = sst.Link("krtr_fwd_port1_" + str(nodeId))
        fwdport1.connect( (rtrfwd, "local1", mesh_link_latency), l3fwd)
        dataport1 = sst.Link("kRtr_data_port1_" + str(nodeId))
        dataport1.connect( (rtrdata, "local1", mesh_link_latency), l3data)


print "Building model..."

kRtrReq=[]
kRtrAck=[]
kRtrFwd=[]
kRtrData=[]
for x in range (0, mesh_stops_x):
    for y in range (0, mesh_stops_y):
        nodeNum = len(kRtrReq)
        kRtrReq.append(sst.Component("krtr_req_" + str(nodeNum), "kingsley.noc_mesh"))
        kRtrReq[-1].addParams(ctrl_network_params)
        kRtrAck.append(sst.Component("krtr_ack_" + str(nodeNum), "kingsley.noc_mesh"))
        kRtrAck[-1].addParams(ctrl_network_params)
        kRtrFwd.append(sst.Component("krtr_fwd_" + str(nodeNum), "kingsley.noc_mesh"))
        kRtrFwd[-1].addParams(ctrl_network_params)
        kRtrData.append(sst.Component("krtr_data_" + str(nodeNum), "kingsley.noc_mesh"))
        kRtrData[-1].addParams(data_network_params)
        
        kRtrReq[-1].addParams({"local_ports" : 2})
        kRtrAck[-1].addParams({"local_ports" : 2})
        kRtrFwd[-1].addParams({"local_ports" : 2})
        kRtrData[-1].addParams({"local_ports" : 2})

if not quiet:
    print "Routers built, connecting..."
i = 0
for y in range(0, mesh_stops_y):
    for x in range (0, mesh_stops_x):

        # North-south connections
        if y != (mesh_stops_y -1):
            kRtrReqNS = sst.Link("krtr_req_ns_" + str(i))
            kRtrReqNS.connect( (kRtrReq[i], "south", mesh_link_latency), (kRtrReq[i + mesh_stops_x], "north", mesh_link_latency) )
            kRtrAckNS = sst.Link("krtr_ack_ns_" + str(i))
            kRtrAckNS.connect( (kRtrAck[i], "south", mesh_link_latency), (kRtrAck[i + mesh_stops_x], "north", mesh_link_latency) )
            kRtrFwdNS = sst.Link("krtr_fwd_ns_" + str(i))
            kRtrFwdNS.connect( (kRtrFwd[i], "south", mesh_link_latency), (kRtrFwd[i + mesh_stops_x], "north", mesh_link_latency) )
            kRtrDataNS = sst.Link("krtr_data_ns_" + str(i))
            kRtrDataNS.connect( (kRtrData[i], "south", mesh_link_latency), (kRtrData[i + mesh_stops_x], "north", mesh_link_latency) )

        if x != (mesh_stops_x - 1):
            kRtrReqEW = sst.Link("krtr_req_ew_" + str(i))
            kRtrReqEW.connect( (kRtrReq[i], "east", mesh_link_latency), (kRtrReq[i+1], "west", mesh_link_latency) )
            kRtrAckEW = sst.Link("krtr_ack_ew_" + str(i))
            kRtrAckEW.connect( (kRtrAck[i], "east", mesh_link_latency), (kRtrAck[i+1], "west", mesh_link_latency) )
            kRtrFwdEW = sst.Link("krtr_fwd_ew_" + str(i))
            kRtrFwdEW.connect( (kRtrFwd[i], "east", mesh_link_latency), (kRtrFwd[i+1], "west", mesh_link_latency) )
            kRtrDataEW = sst.Link("krtr_data_ew_" + str(i))
            kRtrDataEW.connect( (kRtrData[i], "east", mesh_link_latency), (kRtrData[i+1], "west", mesh_link_latency) )

        setNode(i, kRtrReq[i], kRtrAck[i], kRtrFwd[i], kRtrData[i])
        i = i + 1

# Enable SST Statistics Outputs for this simulation
sst.setStatisticLoadLevel(16)
sst.enableAllStatisticsForAllComponents({"type":"sst.AccumulatorStatistic"})

sst.setStatisticOutput("sst.statOutputCSV")
sst.setStatisticOutputOptions( {
        "filepath"  : statFile,
        "separator" : ", "
} )
print "Model complete"
