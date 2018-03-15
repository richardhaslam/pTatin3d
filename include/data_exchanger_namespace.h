
/*
 Set of macros to namespace DataExchanger symbols and objects.
 This is required for compatibility with PETSc v3.8
*/

#ifndef __PTATIN_DATA_EXCHANGER_NAMESPACE_H__
#define __PTATIN_DATA_EXCHANGER_NAMESPACE_H__

/* contants */
#define status_names ptatin_status_names

/* enums */
#define DEObjectState ptatinDEObjectState

/* types */
#define _p_DataEx _p_ptatinDataEx
#define DataEx    ptatinDataEx

/* symbols */

#define DataExCreate                ptatinDataExCreate
#define DataExView                  ptatinDataExView
#define DataExDestroy               ptatinDataExDestroy
#define DataExTopologyInitialize    ptatinDataExTopologyInitialize
#define DataExTopologyAddNeighbour  ptatinDataExTopologyAddNeighbour
#define DataExTopologyFinalize      ptatinDataExTopologyFinalize
#define DataExInitializeSendCount   ptatinDataExInitializeSendCount
#define DataExAddToSendCount        ptatinDataExAddToSendCount
#define DataExFinalizeSendCount     ptatinDataExFinalizeSendCount
#define DataExPackInitialize        ptatinDataExPackInitialize
#define DataExPackData              ptatinDataExPackData
#define DataExPackFinalize          ptatinDataExPackFinalize
#define DataExBegin                 ptatinDataExBegin
#define DataExEnd                   ptatinDataExEnd
#define DataExGetSendData           ptatinDataExGetSendData
#define DataExGetRecvData           ptatinDataExGetRecvData
#define DataExTopologyGetNeighbours ptatinDataExTopologyGetNeighbours

#endif
