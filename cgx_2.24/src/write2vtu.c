/* --------------------------------------------------------------------  */
/* write2vtu.c - Convert CalculiX results to VTK XML (VTU) Binary        */
/* MODIFIED VERSION:                                                     */
/* - RAW OUTPUT (No scaling, No unit conversion)                         */
/* - FIXED MESH TOPOLOGY (Node ID Mapping included)                      */
/* - FIXED PE15 ELEMENT TYPE (15-node Wedge)                             */
/* --------------------------------------------------------------------  */

#include <extUtil.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
void calcDatasets( int lc, Summen *anz, Nodes *node, Datasets *lcase );

/* Helper: number of nodes per element type */
int get_elem_node_count(int type) {
    switch (type) {
        case 1:  return 8;  // Hex8
        case 2:  return 6;  // Wedge6
        case 3:  return 4;  // Tet4
        case 4:  return 20; // Hex20
        case 5:  return 15; // pe15 -> Quadratic Wedge (15 nodes)
        case 6:  return 10; // Tet10
        case 7:  return 3;  // Tri3
        case 8:  return 6;  // Tri6
        case 9:  return 4;  // Quad4
        case 10: return 8;  // Quad8
        case 11: return 2;  // Line2
        case 12: return 3;  // Line3
        default: return 0;
    }
}

/* Helper: CalculiX type -> VTK cell type */
int get_vtk_cell_type(int type) {
    switch (type) {
        case 1:  return 12; // VTK_HEXAHEDRON
        case 2:  return 13; // VTK_WEDGE
        case 3:  return 10; // VTK_TETRA
        case 4:  return 25; // VTK_QUADRATIC_HEXAHEDRON
        case 5:  return 26; // VTK_QUADRATIC_WEDGE (Corrected from Pyramid)
        case 6:  return 24; // VTK_QUADRATIC_TETRA
        case 7:  return 5;  // VTK_TRIANGLE
        case 8:  return 22; // VTK_QUADRATIC_TRIANGLE
        case 9:  return 9;  // VTK_QUAD
        case 10: return 23; // VTK_QUADRATIC_QUAD
        case 11: return 3;  // VTK_LINE
        case 12: return 21; // VTK_QUADRATIC_EDGE
        default: return 0;
    }
}

int write2vtu(char *datout, Summen *anz, Nodes *node, Elements *elem, Datasets *lcase, char **dat)
{
    FILE *handle = NULL;
    int i, j, lc;
    int node_count = anz->n; // Actual number of nodes
    int elem_count = anz->e;
    
    // --- MAPOWANIE W?Z?ÓW ---
    int *nodeMap = (int*)malloc((anz->nmax + 1) * sizeof(int));
    if (nodeMap == NULL) {
        printf("ERROR: Out of memory for node map.\n");
        return -1;
    }
    for(i = 0; i <= anz->nmax; i++) nodeMap[i] = -1;
    
    // Wype?nienie mapy
    for(i = 0; i < anz->n; i++) {
        if(node[i].nr <= anz->nmax) {
            nodeMap[node[i].nr] = i;
        }
    }

    /* Offset counters for Appended mode */
    uint64_t current_offset = 0;
    uint64_t offset_points = 0;
    uint64_t offset_connectivity = 0;
    uint64_t offset_offsets = 0;
    uint64_t offset_types = 0;

    /* Block sizes */
    uint64_t size_points_block = (uint64_t)node_count * 3 * sizeof(double);
    uint64_t size_conn_block = 0;
    uint64_t size_offsets_block = (uint64_t)elem_count * sizeof(int32_t);
    uint64_t size_types_block = (uint64_t)elem_count * sizeof(uint8_t);

    /* Calculate connectivity size */
    for (i = 0; i < elem_count; i++) {
        size_conn_block += get_elem_node_count(elem[i].type) * sizeof(int32_t);
    }

    printf("Exporting to VTU (Binary/Raw). No scaling applied.\n");

    /* Safe output filename */
    char *outfile = malloc(strlen(datout) + 5);
    if (!outfile) { free(nodeMap); return -1; }
    sprintf(outfile, "%s.vtu", datout);

    handle = fopen(outfile, "wb");
    free(outfile);
    if (handle == NULL) {
        printf("\nThe output file \"%s.vtu\" could not be opened.\n\n", datout);
        free(nodeMap);
        return -1;
    }
    printf(" file %s.vtu opened\n", datout);

    /* --- XML Header --- */
    fprintf(handle, "<?xml version=\"1.0\"?>\n");
    fprintf(handle, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
    fprintf(handle, "  <UnstructuredGrid>\n");
    fprintf(handle, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", node_count, elem_count);

    /* --- FieldData: Time values --- */
    fprintf(handle, "      <FieldData>\n");
    for (lc = 0; lc < anz->l; lc++) {
        fprintf(handle, "        <DataArray type=\"Float64\" Name=\"TimeValue_%d\" NumberOfComponents=\"1\" format=\"ascii\">\n", lc);
        fprintf(handle, "          %g\n", lcase[lc].value);
        fprintf(handle, "        </DataArray>\n");
    }
    fprintf(handle, "      </FieldData>\n");

    /* --- Points --- */
    offset_points = current_offset;
    current_offset += sizeof(uint32_t) + size_points_block;

    fprintf(handle, "      <Points>\n");
    fprintf(handle, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lu\"/>\n",
            (unsigned long)offset_points);
    fprintf(handle, "      </Points>\n");

    /* --- Cells --- */
    offset_connectivity = current_offset;
    current_offset += sizeof(uint32_t) + size_conn_block;

    offset_offsets = current_offset;
    current_offset += sizeof(uint32_t) + size_offsets_block;

    offset_types = current_offset;
    current_offset += sizeof(uint32_t) + size_types_block;

    fprintf(handle, "      <Cells>\n");
    fprintf(handle, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%lu\"/>\n",
            (unsigned long)offset_connectivity);
    fprintf(handle, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%lu\"/>\n",
            (unsigned long)offset_offsets);
    fprintf(handle, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"%lu\"/>\n",
            (unsigned long)offset_types);
    fprintf(handle, "      </Cells>\n");

    /* --- PointData --- */
    fprintf(handle, "      <PointData>\n");

    for (lc = 0; lc < anz->l; lc++) {
        if (!lcase[lc].loaded) {
            if (readfrdblock(lc, anz, node, lcase) == -1) {
                printf("ERROR: Could not read dataset %d\n", lc + 1);
                continue;
            }
            calcDatasets(lc, anz, node, lcase);
        }

        char base_name[64] = "";
        int comps = 0;
        int is_tensor = 0;
        int is_contact = 0;
        // is_temp usuni?te, bo nie skalujemy ju? temperatury

        if (compare(lcase[lc].name, "DISP", 4) == 4) {
            strcpy(base_name, "DISP"); comps = 3;
        } else if (compare(lcase[lc].name, "FORC", 4) == 4) {
            strcpy(base_name, "FORC"); comps = 3;
        } else if (compare(lcase[lc].name, "TEMP", 4) == 4 || compare(lcase[lc].name, "NDTEMP", 6) == 6) {
            strcpy(base_name, "TEMP"); comps = 1; 
        } else if (compare(lcase[lc].name, "STRESS", 6) == 6) {
            strcpy(base_name, "STRESS"); comps = 6; is_tensor = 1;
        } else if (compare(lcase[lc].name, "TOSTRAIN", 8) == 8 || compare(lcase[lc].name, "STRAIN", 6) == 6) {
            strcpy(base_name, "STRAIN"); comps = 6; is_tensor = 1;
        } else if (compare(lcase[lc].name, "CONTACT", 7) == 7) {
            strcpy(base_name, "CPRESS"); comps = 1;
        } else if (compare(lcase[lc].name, "PE", 2) == 2) {
            strcpy(base_name, "PE"); comps = 1;
        }

        if (comps > 0) {
            char comp_names[1024] = "";
            if (comps == 6) {
                snprintf(comp_names, sizeof(comp_names),
                    " ComponentName0=\"XX\" ComponentName1=\"YY\" ComponentName2=\"ZZ\""
                    " ComponentName3=\"XY\" ComponentName4=\"YZ\" ComponentName5=\"ZX\"");
            }

            uint64_t ds_size = (uint64_t)node_count * comps * sizeof(double);
            fprintf(handle, "        <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"appended\" offset=\"%lu\"%s/>\n",
                    base_name, comps, (unsigned long)current_offset, comp_names);
            current_offset += sizeof(uint32_t) + ds_size;

            if (is_tensor) {
                fprintf(handle, "        <DataArray type=\"Float64\" Name=\"vMises_%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lu\"/>\n",
                        base_name, (unsigned long)current_offset);
                current_offset += sizeof(uint32_t) + (uint64_t)node_count * sizeof(double);

                fprintf(handle, "        <DataArray type=\"Float64\" Name=\"sgMises_%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lu\"/>\n",
                        base_name, (unsigned long)current_offset);
                current_offset += sizeof(uint32_t) + (uint64_t)node_count * sizeof(double);
            }
        }

        /* SDV */
        if (compare(lcase[lc].name, "SDV", 3) == 3 && lcase[lc].ncomps > 0) {
            for (int c = 0; c < lcase[lc].ncomps; c++) {
                const char *comp_name = (lcase[lc].compName[c] && strlen(lcase[lc].compName[c])) ?
                                        lcase[lc].compName[c] : "UNKNOWN";
                uint64_t ds_size = (uint64_t)node_count * sizeof(double);
                fprintf(handle, "        <DataArray type=\"Float64\" Name=\"%s_%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lu\"/>\n",
                        lcase[lc].name, comp_name, (unsigned long)current_offset);
                current_offset += sizeof(uint32_t) + ds_size;
            }
        }
    }

    fprintf(handle, "      </PointData>\n");
    fprintf(handle, "    </Piece>\n");
    fprintf(handle, "  </UnstructuredGrid>\n");

    /* --- Appended Data --- */
    fprintf(handle, "  <AppendedData encoding=\"raw\">\n");
    fprintf(handle, "_");

    /* 1. POINTS (RAW DATA) */
    {
        uint32_t block_len = (uint32_t)size_points_block;
        fwrite(&block_len, sizeof(uint32_t), 1, handle);
        double *pos = malloc(node_count * 3 * sizeof(double));
        if (!pos) { free(nodeMap); fclose(handle); return -1; }
        
        for (i = 0; i < node_count; i++) {
            // BEZ SKALOWANIA
            pos[3*i]     = node[node[i].nr].nx;
            pos[3*i + 1] = node[node[i].nr].ny;
            pos[3*i + 2] = node[node[i].nr].nz;
        }
        fwrite(pos, sizeof(double), node_count * 3, handle);
        free(pos);
    }

    /* 2. CONNECTIVITY (With Node Mapping) */
    {
        uint32_t block_len = (uint32_t)size_conn_block;
        fwrite(&block_len, sizeof(uint32_t), 1, handle);
        int32_t *conn = malloc(size_conn_block);
        if (!conn) { free(nodeMap); fclose(handle); return -1; }
        
        int idx = 0;
        for (i = 0; i < elem_count; i++) {
            int n = get_elem_node_count(elem[i].type);
            int *ids = elem[i].nod;
            
            if (elem[i].type == 4) { // Hex20 (Permutation needed)
                for (j = 0; j < 12; j++)  conn[idx++] = nodeMap[ids[j]];
                for (j = 16; j < 20; j++) conn[idx++] = nodeMap[ids[j]];
                for (j = 12; j < 16; j++) conn[idx++] = nodeMap[ids[j]];
            } else {
                for (j = 0; j < n; j++) {
                    conn[idx++] = nodeMap[ids[j]];
                }
            }
        }
        fwrite(conn, 1, size_conn_block, handle);
        free(conn);
    }

    /* 3. OFFSETS */
    {
        uint32_t block_len = (uint32_t)size_offsets_block;
        fwrite(&block_len, sizeof(uint32_t), 1, handle);
        int32_t *off = malloc(elem_count * sizeof(int32_t));
        if (!off) { free(nodeMap); fclose(handle); return -1; }
        
        int32_t accum = 0;
        for (i = 0; i < elem_count; i++) {
            accum += get_elem_node_count(elem[i].type);
            off[i] = accum;
        }
        fwrite(off, sizeof(int32_t), elem_count, handle);
        free(off);
    }

    /* 4. TYPES */
    {
        uint32_t block_len = (uint32_t)size_types_block;
        fwrite(&block_len, sizeof(uint32_t), 1, handle);
        uint8_t *types = malloc(elem_count * sizeof(uint8_t));
        if (!types) { free(nodeMap); fclose(handle); return -1; }
        
        for (i = 0; i < elem_count; i++) {
            types[i] = (uint8_t)get_vtk_cell_type(elem[i].type);
        }
        fwrite(types, sizeof(uint8_t), elem_count, handle);
        free(types);
    }

    /* 5. POINT DATA VALUES (RAW) */
    for (lc = 0; lc < anz->l; lc++) {
        if (!lcase[lc].loaded) continue;

        char base_name[64] = "";
        int comps = 0;
        int is_tensor = 0;
        int is_contact = 0;

        if (compare(lcase[lc].name, "DISP", 4) == 4) { strcpy(base_name, "DISP"); comps = 3; }
        else if (compare(lcase[lc].name, "FORC", 4) == 4) { strcpy(base_name, "FORC"); comps = 3; }
        else if (compare(lcase[lc].name, "TEMP", 4) == 4 || compare(lcase[lc].name, "NDTEMP", 6) == 6) { strcpy(base_name, "TEMP"); comps = 1; }
        else if (compare(lcase[lc].name, "STRESS", 6) == 6) { strcpy(base_name, "STRESS"); comps = 6; is_tensor = 1; }
        else if (compare(lcase[lc].name, "TOSTRAIN", 8) == 8 || compare(lcase[lc].name, "STRAIN", 6) == 6) { strcpy(base_name, "STRAIN"); comps = 6; is_tensor = 1; }
        else if (compare(lcase[lc].name, "CONTACT", 7) == 7) { strcpy(base_name, "CPRESS"); comps = 1; is_contact = 1; }
        else if (compare(lcase[lc].name, "PE", 2) == 2) { strcpy(base_name, "PE"); comps = 1; }

        if (comps > 0) {
            uint32_t block_len = (uint32_t)(node_count * comps * sizeof(double));
            fwrite(&block_len, sizeof(uint32_t), 1, handle);

            double *data = malloc(node_count * comps * sizeof(double));
            if (!data) { free(nodeMap); fclose(handle); return -1; }

            for (i = 0; i < node_count; i++) {
                int nid = node[i].nr;
                
                if (comps == 3) {
                    data[3*i]     = lcase[lc].dat[0][nid];
                    data[3*i + 1] = lcase[lc].dat[1][nid];
                    data[3*i + 2] = lcase[lc].dat[2][nid];
                } else if (comps == 1) {
                    // BEZ KONWERSJI TEMP (zapisujemy warto?? bezpo?rednio)
                    double val = is_contact ? lcase[lc].dat[3][nid] : lcase[lc].dat[0][nid];
                    data[i] = val; 
                } else if (comps == 6) {
                    data[6*i + 0] = lcase[lc].dat[0][nid]; 
                    data[6*i + 1] = lcase[lc].dat[1][nid]; 
                    data[6*i + 2] = lcase[lc].dat[2][nid]; 
                    data[6*i + 3] = lcase[lc].dat[3][nid]; 
                    data[6*i + 4] = lcase[lc].dat[5][nid]; 
                    data[6*i + 5] = lcase[lc].dat[4][nid]; 
                }
            }
            fwrite(data, sizeof(double), node_count * comps, handle);
            free(data);

            if (is_tensor) {
                // vMises
                uint32_t scalar_len = (uint32_t)(node_count * sizeof(double));
                fwrite(&scalar_len, sizeof(uint32_t), 1, handle);
                double *vm = malloc(node_count * sizeof(double));
                for (i = 0; i < node_count; i++) vm[i] = lcase[lc].dat[6][node[i].nr];
                fwrite(vm, sizeof(double), node_count, handle);
                free(vm);

                // sgMises
                fwrite(&scalar_len, sizeof(uint32_t), 1, handle);
                double *sg = malloc(node_count * sizeof(double));
                for (i = 0; i < node_count; i++) sg[i] = lcase[lc].dat[25][node[i].nr];
                fwrite(sg, sizeof(double), node_count, handle);
                free(sg);
            }
        }

        /* SDV */
        if (compare(lcase[lc].name, "SDV", 3) == 3 && lcase[lc].ncomps > 0) {
            uint32_t scalar_len = (uint32_t)(node_count * sizeof(double));
            for (int c = 0; c < lcase[lc].ncomps; c++) {
                fwrite(&scalar_len, sizeof(uint32_t), 1, handle);
                double *buf = malloc(node_count * sizeof(double));
                if (!buf) { free(nodeMap); fclose(handle); return -1; }
                for (i = 0; i < node_count; i++) {
                    buf[i] = lcase[lc].dat[c][node[i].nr];
                }
                fwrite(buf, sizeof(double), node_count, handle);
                free(buf);
            }
        }

        if (lcase[lc].loaded) {
            for (int k = 0; k < lcase[lc].ncomps; k++) free(lcase[lc].dat[k]);
            free(lcase[lc].dat);
            lcase[lc].dat = NULL;
            lcase[lc].loaded = 0;
        }
    }

    fprintf(handle, "\n  </AppendedData>\n");
    fprintf(handle, "</VTKFile>\n");

    free(nodeMap);
    fclose(handle);
    return 1;
}
